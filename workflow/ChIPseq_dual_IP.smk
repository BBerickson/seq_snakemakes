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

INDEXES = config['INDEXES']
if isinstance(INDEXES, str):
    INDEXES = [INDEXES]

# Load genome-specific config files into a dictionary
config_indexes = {}
for idx in INDEXES:
    path = Path("workflow/ref") / f"{idx}.yaml"
    if not path.exists():
        sys.exit(f"ERROR: {idx} is not a valid genome index.")
    with open(path) as f:
        config_indexes[idx] = yaml.safe_load(f)


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
USER         = config.get("USER")
UNALIGNED    = config.get("UNALIGNED")

# Directories for data and scripts
FASTQ_DIR = PROJ + "/fastqs"

os.makedirs(FASTQ_DIR, exist_ok = True)

# Simplify ALL_SAMPLES dictionary
SAMPLES, SAMPIN, GROUPS, NORMMAP, PAIREDMAP = process_samples(
    ALL_SAMPLES, INDEXES, NORM, ORIENTATION
)

# align to the second entry then the first
INDEX_SECOND = INDEXES[0]
INDEX_PATH   = config_indexes.get(INDEXES[0], {}).get("INDEX_PATH", config["INDEX_PATH"])
INDEX_PATH2  = config_indexes.get(INDEXES[1], {}).get("INDEX_PATH", config["INDEX_PATH"])
INDEX_MAP    = config_indexes.get(INDEXES[1], {}).get("INDEX_MAP", config["INDEX_MAP"])

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
    NORMMAP[key] = [(index, norm, f'scalefactor_{INDEXES[-1]}' if suffix.lower() == 'scalefactor' else suffix) 
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

# for optinal subsetting bam directory
BAM_PATH = _get_bampath(NORM)
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
        expand("{proj}/stats/{proj}_{step}.tsv", proj=PROJ, step=["fastqc", "clumpify", "bbduk", "aligned"]),
        
        # bam URL
        [] if config.get("skip_bam_urls") else [
            expand(PROJ + "/URLS/" + PROJ + "_{index}_bam_URL.txt", index=INDEXES)
        ],
        
        # results
        [] if config.get("stop_after_alignment") else [
          expand(PROJ + "/stats/" + PROJ + "_{index}_fragment_results.tsv", index=INDEXES),
          expand(PROJ + "/report/" + PROJ + "_{index}_fragmentSize.pdf", index=INDEXES)
        ],
        PROJ + "/report/" + PROJ + "_results.tsv",
        [] if config.get("skip_html_report") else [
            expand(SEQ_DATE + "_" + PROJ + "_{index}_qc_analysis.html", index=INDEXES[0])
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
