# ===== Snake file for making bidirectional bigwig + matrix files ================

# Configure shell for all rules
shell.executable("/bin/bash")
shell.prefix("[ -f ~/.bash_profile ] && source ~/.bash_profile; set -o nounset -o pipefail -o errexit -x; ")

# Python packages
import os
import sys
import yaml
from pathlib import Path

# Include custom Python functions
include: workflow.source_path("scripts/funs.py")

# ------------------------------------------------------------------------------
# Load index config
# ------------------------------------------------------------------------------

INDEXES, config_indexes, ROLE_MAP, ALIGNMENT_MODE, PRIMARY_INDEX, INDEX_PATH, INDEX_MAP = (
    load_index_config(config)
)
# MASK is a bamCoverage/bigwig-level concern -- bigwigs don't vary by matrix
# region subset, so MASK always comes from the primary INDEXES entry's own
# yaml, never from a GENELISTS variant below.
MASK        = config_indexes[PRIMARY_INDEX].get("MASK")
SPIKE_INDEX = ROLE_MAP.get("spike")   # None for single mode
INDEXES     = [PRIMARY_INDEX]         # matrix pipelines stay architecturally single-genome downstream

# GENELISTS: one or more gene-list variants of the primary index to build
# matrix files for in this one run (e.g. whole-genome + tRNA-only), keyed by
# each variant's own GENELIST tag. Defaults to just the primary index's own
# yaml when samples.yaml has no GENELISTS -- identical to today's
# single-genelist behavior.
GENELIST_VARIANTS = load_genelists(config, config_indexes, PRIMARY_INDEX)

# ------------------------------------------------------------------------------
# Assign parameters from configs
# ------------------------------------------------------------------------------

# Docker container
SINGULARITY = config.get("CONTAINER", "")
# Check if we should use Singularity
USE_SINGULARITY = bool(SINGULARITY and SINGULARITY.strip())
# Set container if using Singularity
if USE_SINGULARITY:
    singularity:
        config["CONTAINER"]

# start_from_bw: fetch/symlink existing bigwigs instead of generating them from
# BAMs via bamCoverage.
START_FROM_BW = _coerce_bool(config.get("START_FROM_BW", False))

# From main config
PROJ         = config.get("PROJ")
RAW_DATA     = config.get("RAW_DATA")
ALL_SAMPLES  = config.get("SAMPLES")
SEQ_DATE     = config.get("SEQ_DATE")
CHROM_SIZES  = config.get("CHROM_SIZES")
CMD_PARAMS   = config.get("CMD_PARAMS")
COLORS       = config.get("COLORS")
NORM         = config.get("NORM")
ORIENTATION  = config.get("ORIENTATION")
PAIRED       = config.get("PAIRED")
REGIONS      = config.get("REGIONS")
USER         = config.get("USER")

if START_FROM_BW:
    BW_DIR = PROJ + "/raw_bw"
    os.makedirs(BW_DIR, exist_ok = True)

# Simplify ALL_SAMPLES dictionary
SAMPLES, SAMPIN, GROUPS, NORMMAP, PAIREDMAP = process_samples(
    ALL_SAMPLES, PRIMARY_INDEX, [], NORM, ORIENTATION, paired=PAIRED
)

# make file suffix from bamCoverage settings and NORM
for sample, norm_list in NORMMAP.items():
    updated_list = [
        (index, norm_value, _get_normtype(
            CMD_PARAMS["bamCoverage"],
            norm_value,
            CMD_PARAMS.get("bamCoverageBL", False),
            ORIENTATION
        ))
        for index, norm_value in norm_list
    ]
    NORMMAP[sample] = updated_list

# add scalefactor index info
NORMMAP = tag_scalefactor_suffix(NORMMAP, SPIKE_INDEX)

# unpack samples and groups
SAMS = [[y, x] for y in SAMPLES for x in SAMPLES[y]]
NAMS = [x[0] for x in SAMS] # newnames
SAMS = [x[1] for x in SAMS] # samples
GRPS = [[y, x] for y in GROUPS for x in GROUPS[y]]
GRPS = [x[0] for x in GRPS] # groups
NAMS_UNIQ = list(dict.fromkeys(NAMS))
GRPS_UNIQ = list(dict.fromkeys(GRPS))
SAMS_UNIQ = list(dict.fromkeys(SAMS))

# Print summary of samples and groups
print("NORMMAP (%s): %s\n" % (len(NORMMAP), NORMMAP))

# Wildcard constraints
WILDCARD_REGEX = r"[a-zA-Z0-9_\-]+" # Matches alphanumeric characters, underscores, and hyphens

wildcard_constraints:
    sample = WILDCARD_REGEX,
    newnam = WILDCARD_REGEX,
    group  = WILDCARD_REGEX,
    index  = WILDCARD_REGEX,
    suffix = WILDCARD_REGEX,
    covarg = r"[a-zA-Z0-9_.-]+",
    region = "543|5|3|PI|EI"

HEATMAP_REGIONS = ["543","5","5L","3"]

COLS_DICT = _get_colors(NAMS_UNIQ, COLORS)

# Build region x covarg combos for every gene-list variant of the primary
# index, and a (region, covarg) -> yaml_data lookup so the matrix rules can
# find the right FW_REF/REV_REF/FW_PI_REF/REV_PI_REF for whichever variant
# produced a given covarg (each variant's GENELIST tag keeps its covarg
# strings distinct -- see load_genelists's uniqueness check).
REF_BY_REGION_COVARG = {}
product = []
for genelist_tag, yaml_data in GENELIST_VARIANTS.items():
    covargs = _get_all_matrixtypes(REGIONS, CMD_PARAMS, genelist_tag)
    for region, covarg in zip(REGIONS, covargs):
        REF_BY_REGION_COVARG[(region, covarg)] = yaml_data
    product += [(s, i, v) for s in NAMS_UNIQ for i, v in zip(REGIONS, covargs)]

# Convert to DataFrame
REGIONS_COVARGS = pd.DataFrame(product, columns=['Newnam', 'Region', 'Value'])
# Combine into a list of records with expanded NormMap
SAM_NORM = []
for key in SAMPIN:
    sam_value = SAMPIN[key][0]
    for index, norm, suffix in NORMMAP[key]:
        bam_dir = "bams_sub" if norm == "subsample" else "bams"
        SAM_NORM.append([sam_value, key, index, norm, suffix, bam_dir])

# Create DataFrame
DF_SAM_NORM = pd.DataFrame(SAM_NORM, columns=['Sample', 'Newnam', 'Index', 'Norm', 'Suffix', 'Bam_dir'])
# Merge on 'Newnam' and 'Index'
DF_SAM_NORM = REGIONS_COVARGS.merge(DF_SAM_NORM, on=['Newnam'], how='left')
# PI matix files don't get normalized
mask = DF_SAM_NORM["Region"] == "PI"
DF_SAM_NORM.loc[mask, "Suffix"] = [
    suffix.replace(norm, "none", 1)
    for suffix, norm in zip(DF_SAM_NORM.loc[mask, "Suffix"], DF_SAM_NORM.loc[mask, "Norm"])
]
DF_SAM_NORM.loc[mask, "Norm"] = "none"

# Keep SAM_NORM in sync with the PI-relabeled DF_SAM_NORM above -- 04a_bamCoverage_stranded.snake's
# embedded scale-factor lookup reads params.sam_norm (this list), and needs the
# PI-region "none"-relabeled rows to exist here too, or it can't find a matching
# row for a PI-region bigwig's suffix and crashes.
SAM_NORM = (
    DF_SAM_NORM[['Sample', 'Newnam', 'Index', 'Norm', 'Suffix', 'Bam_dir']]
    .drop_duplicates()
    .values.tolist()
)

print(DF_SAM_NORM.to_string(index=False))

# Final output files
rule all:
    input:
        # bamCoverage
        [] if config.get("stop_after_alignment") else [
          expand(
              PROJ + "/bw/{newnam}_aligned_{index}_" + SEQ_DATE + "_norm_{suffix}_fw.bw",
              zip,
              newnam=DF_SAM_NORM['Newnam'],
              index=DF_SAM_NORM['Index'],
              suffix=DF_SAM_NORM['Suffix']
          ),
          expand(
              PROJ + "/bw/{newnam}_aligned_{index}_" + SEQ_DATE + "_norm_{suffix}_rev.bw",
              zip,
              newnam=DF_SAM_NORM['Newnam'],
              index=DF_SAM_NORM['Index'],
              suffix=DF_SAM_NORM['Suffix']
          ),

          # bamCoverage URL for amc-sandbox (only relevant when generating
          # bigwigs here; when they come from an already-completed alignment
          # pipeline, its own run already made these URLs)
          [] if (not START_FROM_BW) or config.get("skip_bw_urls") else [
              expand(
                  PROJ + "/URLS/" + PROJ + "_{index}_" + SEQ_DATE + "_norm_{suffix}_bw_URL.txt",
                  zip, index=DF_SAM_NORM['Index'], suffix=DF_SAM_NORM['Suffix']
              )
          ],

          # matrix url file for amc-sandbox
          [] if config.get("skip_matrix_url") else [
            expand(
                PROJ + "/URLS/{region}_aligned_{index}_" + SEQ_DATE + "_{covarg}_norm_{suffix}_bidirectonal_matrix.url.txt",
                zip,
                region=DF_SAM_NORM['Region'],
                index=DF_SAM_NORM['Index'],
                covarg=DF_SAM_NORM['Value'],
                suffix=DF_SAM_NORM['Suffix']
            )
          ]
        ]

# bigwigs: generate from BAMs, or fetch existing ones
if START_FROM_BW:
    include: "rules/04_get_BW_Stranded.snake"
    include: "rules/04b_bw_UCSC_URL_stranded.snake"
else:
    include: "rules/04a_bamCoverage_stranded.snake"
# 5 sense matrix file
include: "rules/05_bidirectional_matrix.snake"
