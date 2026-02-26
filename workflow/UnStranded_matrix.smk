# ===== Snake file for processing ChIP-seq data ================================

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

# path of samples file
SAMPLES_FILE = workflow.configfiles[0] if workflow.configfiles else "samples.yaml"

# ------------------------------------------------------------------------------
# Load main genome config
# ------------------------------------------------------------------------------

GENOME = config["GENOME"]
GENOME_CONFIG = Path("workflow/ref") / f"{GENOME}.yaml"

if not GENOME_CONFIG.exists():
    sys.exit(f"ERROR: {GENOME} is not a valid GENOME selection.")

# Load the main config file into Snakemake's config dictionary
configfile: str(GENOME_CONFIG)

# ------------------------------------------------------------------------------
# Load additional genome-specific configs manually
# ------------------------------------------------------------------------------

raw_indexes = config['INDEXES']
INDEXES = [raw_indexes] if isinstance(raw_indexes, str) else [raw_indexes[0]]
INDEXES_LAST = [raw_indexes] if isinstance(raw_indexes, str) else [raw_indexes[-1]]

# Paths to additional config files
if config.get("GENELIST"):
    GENOME_CONFIG1 = GENOME_CONFIG
else:
    GENOME_CONFIG1 = Path("workflow/ref") / f"{INDEXES[0]}.yaml"

# Validate existence
if not GENOME_CONFIG1.exists():
    sys.exit(f"ERROR: Config file not found for index '{INDEXES[0]}'. Expected at: {GENOME_CONFIG1}")

# Load additional configs
with open(GENOME_CONFIG1) as f:
    config1 = yaml.safe_load(f)

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
INDEX_PATH   = config.get("INDEX_PATH")
INDEX_MAP    = config.get("INDEX_MAP")
NORM         = config.get("NORM")
CMD_PARAMS   = config.get("CMD_PARAMS")
COLORS       = config.get("COLORS")
ORIENTATION  = config.get("ORIENTATION")
SENSE_ASENSE = config.get("SENSE_ASENSE", [])
REGIONS      = config.get("REGIONS")
USER         = config.get("USER")

if not SENSE_ASENSE:
    SENSE_ASENSE = [""]  # Empty string for no sense/antisense distinction
    
# From additional configs
MY_REF      = config1.get("MY_REF")
PI_REF      = config1.get("PI_REF")
GENELIST = config1.get("GENELIST") or ""


# Simplify ALL_SAMPLES dictionary
SAMPLES, SAMPIN, GROUPS, NORMMAP, PAIREDMAP = process_samples(
    ALL_SAMPLES, INDEXES, NORM, ORIENTATION
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
for key in NORMMAP:
    NORMMAP[key] = [(index, norm, f'scalefactor_{INDEXES_LAST[0]}' if suffix.lower() == 'scalefactor' else suffix) 
                    for index, norm, suffix in NORMMAP[key]]

# Combine into a list of records with expanded NormMap
SAM_NORM = []
for key in SAMPIN:
    sam_value = SAMPIN[key][0]
    for index, norm, suffix in NORMMAP[key]:
        SAM_NORM.append([sam_value, key, index, norm, suffix])

# Create DataFrame
DF_SAM_NORM = pd.DataFrame(SAM_NORM, columns=['Sample', 'Newnam', 'Index', 'Norm', 'Suffix'])

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
print(DF_SAM_NORM.to_string(index=False))

# Wildcard constraints
WILDCARD_REGEX = r"[a-zA-Z0-9_\-]+" # Matches alphanumeric characters, underscores, and hyphens

wildcard_constraints:
    sample = WILDCARD_REGEX,
    newnam = WILDCARD_REGEX,
    group  = WILDCARD_REGEX,
    index  = WILDCARD_REGEX,
    suffix = WILDCARD_REGEX,
    covarg = r"[a-zA-Z0-9_.-]+",
    region = "543|5|5L|3|PI|EI",
    sense_asense = r"(sense|anti|)"  # Allow empty, "sense", or "anti"

HEATMAP_REGIONS = ["543","5","5L","3"]

COLS_DICT = _get_colors(NAMS_UNIQ, COLORS)

NORMS = _get_normtype(CMD_PARAMS["bamCoverage"],NORM,CMD_PARAMS.get("bamCoverageBL", False),ORIENTATION)

BAM_PATH = _get_bampath(NORM)

COVARGS = _get_all_matrixtypes(REGIONS,CMD_PARAMS,GENELIST)
GENELIST = config1.get("GENELIST") or "all"

# Create the Cartesian product
product = [(s, i, v) for s in NAMS_UNIQ for i, v in zip(REGIONS, COVARGS)]

# Convert to DataFrame
REGIONS_COVARGS = pd.DataFrame(product, columns=['Newnam', 'Region', 'Value'])
# Combine into a list of records with expanded NormMap
SAM_NORM = []
for key in SAMPIN:
    sam_value = SAMPIN[key][0]
    for index, norm, suffix in NORMMAP[key]:
        SAM_NORM.append([sam_value, key, index, norm, suffix])

# Create DataFrame
DF_SAM_NORM = pd.DataFrame(SAM_NORM, columns=['Sample', 'Newnam', 'Index', 'Norm', 'Suffix'])
print(DF_SAM_NORM.to_string(index=False))

# Merge on 'Newnam' and 'Index'
DF_SAM_NORM = REGIONS_COVARGS.merge(DF_SAM_NORM, on=['Newnam'], how='left')

print(DF_SAM_NORM.to_string(index=False))

# Final output files
rule all:
    input:
        # bamCoverage
        [] if config.get("stop_after_alignment") else [
          expand(
              PROJ + "/bw/{newnam}_aligned_{index}_" + SEQ_DATE + "_norm_{suffix}.bw",
              zip,
              newnam=DF_SAM_NORM['Newnam'],
              index=DF_SAM_NORM['Index'],
              suffix=DF_SAM_NORM['Suffix']
          ),
          
          # matrix files
          expand(
              PROJ + "/matrix/{region}/{newnam}_aligned_{index}_" + SEQ_DATE + "_{region}_{covarg}_norm_{suffix}_matrix.gz",
              zip, 
              region=DF_SAM_NORM['Region'], 
              newnam=DF_SAM_NORM['Newnam'],
              index=DF_SAM_NORM['Index'], 
              covarg=DF_SAM_NORM['Value'], 
              suffix=DF_SAM_NORM['Suffix']
          ),
          
          # matrix url file for amc-sandbox
          [] if config.get("skip_matrix_url") else [
            expand(
                PROJ + "/URLS/{region}_aligned_{index}_" + SEQ_DATE + "_{covarg}_norm_{suffix}_matrix.url.txt",
                zip, 
                region=DF_SAM_NORM['Region'], 
                index=DF_SAM_NORM['Index'], 
                covarg=DF_SAM_NORM['Value'], 
                suffix=DF_SAM_NORM['Suffix']
            )
          ],
          
          # qc with heatmap, cluster, profile plots
          [] if config.get("skip_matrix_html_report") else [
              expand(SEQ_DATE + "_" + PROJ + "_{index}_qc_plots_analysis.html", index=INDEXES[0])
          ]
        ]
        
        
# BW with deeptools bamCoverage
include: "rules/04a_bamCoverage.snake"
# make matrix files
include: "rules/05_UnStranded_matrix.snake"
include: "rules/06_matrix_heatmap.snake"
include: "rules/07_results_html.snake"
