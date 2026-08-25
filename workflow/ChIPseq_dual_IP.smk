# ===== Snake file for processing Bowtie ================================

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
MASK = config_indexes[PRIMARY_INDEX].get("MASK")
INDEX_ROLE = {v: k for k, v in ROLE_MAP.items()}  # index name -> role, for rules/03a_fragmentSize.snake

# sequential mode: two real, separate indexes aligned in two passes. The
# spike-role index is aligned FIRST (bowtie_first_pass); the primary-role
# index is aligned SECOND and is the reported/results index.
SPIKE_INDEX  = ROLE_MAP["spike"]
INDEX_MAP    = config_indexes[SPIKE_INDEX]["INDEX_MAP"]
INDEX_PATH   = config_indexes[SPIKE_INDEX]["INDEX_PATH"]
INDEX_SECOND = config_indexes[PRIMARY_INDEX]["INDEX_MAP"]
INDEX_PATH2  = config_indexes[PRIMARY_INDEX]["INDEX_PATH"]

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

# From main config
PROJ         = config.get("PROJ")
RAW_DATA     = config.get("RAW_DATA")
ALL_SAMPLES  = config.get("SAMPLES")
SEQ_DATE     = config.get("SEQ_DATE")
BARCODES     = config.get("BARCODES")
CMD_PARAMS   = config.get("CMD_PARAMS")
NORM         = config.get('NORM')
COLORS       = config.get("COLORS")
ORIENTATION  = config.get("ORIENTATION")
PAIRED       = config.get("PAIRED")
USER         = config.get("USER")
UNALIGNED    = config.get("UNALIGNED")

if ALIGNMENT_MODE != "sequential":
    sys.exit(f"ERROR: ChIPseq_dual_IP.smk requires ALIGNMENT_MODE: sequential, got {ALIGNMENT_MODE!r}.")

# Directories for data and scripts
FASTQ_DIR = PROJ + "/fastqs"

os.makedirs(FASTQ_DIR, exist_ok = True)

# Simplify ALL_SAMPLES dictionary
SAMPLES, SAMPIN, GROUPS, NORMMAP, PAIREDMAP = process_samples(
    ALL_SAMPLES, PRIMARY_INDEX, [i for i in INDEXES if i != PRIMARY_INDEX], NORM, ORIENTATION, paired=PAIRED
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

# Combine into a list of records with expanded NormMap
SAM_NORM = []
for key in SAMPIN:
    sam_value = SAMPIN[key][0]
    for index, norm, suffix in NORMMAP[key]:
        bam_dir = "bams_sub" if norm == "subsample" else "bams"
        SAM_NORM.append([sam_value, key, index, norm, suffix, bam_dir])

# Create DataFrame
DF_SAM_NORM = pd.DataFrame(SAM_NORM, columns=['Sample', 'Newnam', 'Index', 'Norm', 'Suffix', 'Bam_dir'])

# unpack samples and groups
SAMS = [[y, x] for y in SAMPIN for x in SAMPIN[y]]
NAMS = [x[0] for x in SAMS] # newnames
SAMS = [x[1] for x in SAMS] # samples
GRPS = [[y, x] for y in GROUPS for x in GROUPS[y]]
GRPS = [x[0] for x in GRPS] # groups
NAMS_UNIQ = list(dict.fromkeys(NAMS))
GRPS_UNIQ = list(dict.fromkeys(GRPS))
SAMS_UNIQ = list(dict.fromkeys(SAMS))

# Print summary of samples and groups
print(DF_SAM_NORM.to_string(index=False))

# Create symlinks for fastqs
FASTQS = [_get_fqs(x, RAW_DATA, FASTQ_DIR, paired=PAIREDMAP[x]) for x in SAMS_UNIQ]
FASTQS = sum(FASTQS, [])

# color dictionary 
COLS_DICT = _get_colors(NAMS_UNIQ, COLORS)

ALIGNER = "bowtie2"

# Wildcard constraints
WILDCARD_REGEX = r"[a-zA-Z0-9_\-]+" # Matches alphanumeric characters, underscores, and hyphens

wildcard_constraints:
    sample = WILDCARD_REGEX,
    newnam = WILDCARD_REGEX,
    group  = WILDCARD_REGEX,
    index  = WILDCARD_REGEX,
    suffix = WILDCARD_REGEX


# Final output files
rule all:
    input:
        # process summary files
        expand("{proj}/stats/{proj}_{step}.tsv", proj=PROJ, step=["fastqc", "clumpify", "bbduk"]),
        PROJ + "/stats/" + PROJ + "_" + INDEX_MAP + "_aligned.tsv",
        
        # bam URL
        [] if config.get("skip_bam_urls") else [
            expand(PROJ + "/URLS/" + PROJ + "_{index}_bam_URL.txt", index=INDEXES)
        ],
        
        # results
        [] if config.get("stop_after_alignment") else [
          expand(PROJ + "/stats/" + PROJ + "_{index}_fragment_results.tsv", index=INDEXES),
          expand(PROJ + "/report/" + PROJ + "_{index}_fragmentSize.pdf", index=INDEXES)
        ],
        expand(PROJ + "/report/" + PROJ + "_{index}_results.tsv", index=PRIMARY_INDEX),
        [] if config.get("skip_html_report") else [
            expand(SEQ_DATE + "_" + PROJ + "_{index}" + "_qc_analysis.html", index=PRIMARY_INDEX)
        ],
        
        # bamCoverage
        [] if config.get("stop_after_alignment") else [
          expand(
              PROJ + "/bw/{newnam}_aligned_{index}_" + SEQ_DATE + "_norm_{suffix}.bw",
              zip,
              newnam=DF_SAM_NORM['Newnam'],
              index=DF_SAM_NORM['Index'],
              suffix=DF_SAM_NORM['Suffix']
          ),
          
          # bamCoverage URL for amc-sandbox
          [] if config.get("skip_bw_urls") else [
              expand(
                PROJ + "/URLS/" + PROJ + "_{index}_" + SEQ_DATE + "_norm_{suffix}_bw_URL.txt",
                zip, index=DF_SAM_NORM['Index'], suffix=DF_SAM_NORM['Suffix']
              )
          ]
        ]


# Run FastQC
include: "rules/01_fastqc.snake"
# Run clumpify
include: "rules/01a_clumpify.snake"
# Run bbmerge
include: "rules/01b_bbduk.snake"
# Align reads 
include: "rules/02a_dual_IP_bowtie.snake"
include: "rules/02b_dual_IP_samtools.snake"
include: "rules/02c_align_URLS.snake"
include: "rules/02s_align_subsample.snake"
# Results
include: "rules/03a_featureCounts.snake"
include: "rules/03a_fragmentSize.snake"
include: "rules/03b_results.snake"
# BW with deeptools bamCoverage
include: "rules/04a_bamCoverage.snake"
include: "rules/04b_bw_UCSC_URL.snake"
