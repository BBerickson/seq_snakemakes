# ===== Snake file for processing Ribo-seq (no-UMI variant) ====================

# Configure shell for all rules
shell.executable("/bin/bash")
shell.prefix("[ -f ~/.bash_profile ] && source ~/.bash_profile; set -o nounset -o pipefail -o errexit -x; ")

# python packages
import os
import sys
import yaml
from pathlib import Path

# Include custom Python functions
include: workflow.source_path("scripts/funs.py")

# ------------------------------------------------------------------------------
# Load index config
# ------------------------------------------------------------------------------

# Ribo-seq aligns with STAR -- reuses the same "_star" ref yaml family
# star.smk/RNAseq_samples.yaml already use (e.g. hg38_star.yaml); new optional
# fields (SORTMERNA_DB, BED12, SALMON_INDEX) live on those same ref yamls.
INDEXES, config_indexes, ROLE_MAP, ALIGNMENT_MODE, PRIMARY_INDEX, INDEX_PATH, INDEX_MAP = (
    load_index_config(config, ref_suffix="_star")
)
MASK = config_indexes[PRIMARY_INDEX].get("MASK")
INDEX_ROLE = {v: k for k, v in ROLE_MAP.items()}  # index name -> role, for rules/03a_fragmentSize.snake

# Shared with rules/01b_riboseq_trimgalore.snake -- selects raw fastqs (this
# file) vs. the UMI-extracted fastq (RiboSeq_UMI.smk)
UMI_MODE = False

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
PROJ            = config.get("PROJ")
RAW_DATA        = config.get("RAW_DATA")
ALL_SAMPLES     = config.get("SAMPLES")
SEQ_DATE        = config.get("SEQ_DATE")
CMD_PARAMS      = config.get("CMD_PARAMS")
NORM            = config.get('NORM')
COLORS          = config.get("COLORS")
ORIENTATION     = config.get("ORIENTATION")
PAIRED          = config.get("PAIRED")  # always "False" for Ribo-seq footprints
USER            = config.get("USER")
CONTAMINANT_REF = config.get("CONTAMINANT_REF", "")

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
NAMS = [x[0] for x in SAMS]
SAMS = [x[1] for x in SAMS]
GRPS = [[y, x] for y in GROUPS for x in GROUPS[y]]
GRPS = [x[0] for x in GRPS]
NAMS_UNIQ = list(dict.fromkeys(NAMS))
GRPS_UNIQ = list(dict.fromkeys(GRPS))
SAMS_UNIQ = list(dict.fromkeys(SAMS))

# Print summary of samples and groups
print(DF_SAM_NORM.to_string(index=False))

# Wildcard constraints
WILDCARD_REGEX = r"[a-zA-Z0-9_\-]+"

wildcard_constraints:
    sample = WILDCARD_REGEX,
    newnam = WILDCARD_REGEX,
    group  = WILDCARD_REGEX,
    index  = WILDCARD_REGEX,
    suffix = WILDCARD_REGEX

# Create symlinks for fastqs
FASTQS = [_get_fqs(x, RAW_DATA, FASTQ_DIR, paired=PAIREDMAP[x]) for x in SAMS_UNIQ]
FASTQS = sum(FASTQS, [])

COLS_DICT = _get_colors(NAMS_UNIQ, COLORS)

ALIGNER = "star"

# Final output files
rule all:
    input:
        # process summary files
        expand("{proj}/stats/{proj}_{step}.tsv", proj=PROJ,
               step=["fastqc", "trim_galore", "sortmerna", "strandedness"]),
        PROJ + "/stats/" + PROJ + "_" + INDEX_MAP + "_aligned.tsv",

        # bam URL
        [] if config.get("skip_bam_urls") else [
            expand(PROJ + "/URLS/" + PROJ + "_{index}_bam_URL.txt", index=INDEXES)
        ],

        # results
        [] if config.get("stop_after_alignment") else [
          expand(PROJ + "/stats/" + PROJ + "_{index}_fragment_results.tsv", index=INDEXES),
          expand(PROJ + "/report/" + PROJ + "_{index}_fragmentSize.pdf", index=INDEXES),
          PROJ + "/stats/" + PROJ + "_ribowaltz_psite.tsv",
          PROJ + "/stats/" + PROJ + "_rseqc.tsv",
          PROJ + "/report/" + PROJ + "_" + PRIMARY_INDEX + "_PCA.png",
          PROJ + "/report/" + PROJ + "_" + PRIMARY_INDEX + "_correlation_heatmap.png"
        ],
        expand(PROJ + "/report/" + PROJ + "_{index}_results.tsv", index=PRIMARY_INDEX),
        [] if config.get("skip_html_report") else [
            expand(SEQ_DATE + "_" + PROJ + "_{index}" + "_qc_analysis.html", index=PRIMARY_INDEX)
        ],

        # bamCoverage (single-end footprints -- still stranded, since
        # ribosome footprints reflect transcript strand)
        [] if config.get("stop_after_alignment") else [
          expand(
              PROJ + "/bw/{newnam}_aligned_{index}_" + SEQ_DATE + "_norm_{suffix}_fw.bw",
              zip, newnam=DF_SAM_NORM['Newnam'], index=DF_SAM_NORM['Index'], suffix=DF_SAM_NORM['Suffix']
          ),
          expand(
              PROJ + "/bw/{newnam}_aligned_{index}_" + SEQ_DATE + "_norm_{suffix}_rev.bw",
              zip, newnam=DF_SAM_NORM['Newnam'], index=DF_SAM_NORM['Index'], suffix=DF_SAM_NORM['Suffix']
          ),
          [] if config.get("skip_bw_urls") else [
              expand(
                PROJ + "/URLS/" + PROJ + "_{index}_" + SEQ_DATE + "_norm_{suffix}_bw_URL.txt",
                zip, index=DF_SAM_NORM['Index'], suffix=DF_SAM_NORM['Suffix']
              )
          ]
        ]


# Run FastQC
include: "rules/01_fastqc.snake"
# Trim adapters (Trim Galore)
include: "rules/01b_riboseq_trimgalore.snake"
# Optional contaminant removal (phiX etc.)
include: "rules/01c_riboseq_bbsplit.snake"
# rRNA removal
include: "rules/01d_riboseq_sortmerna.snake"
# Advisory strandedness/library-type inference
include: "rules/01e_riboseq_strandedness.snake"
# Align reads (STAR, genome + transcriptome BAM)
include: "rules/02a_align_star_riboseq.snake"
include: "rules/02b_align_samtools.snake"
include: "rules/02c_align_URLS.snake"
include: "rules/02s_align_mask_subsample.snake"
# Results
include: "rules/03a_featureCounts.snake"
include: "rules/03a_fragmentSize.snake"
include: "rules/03b_riboWaltz.snake"
include: "rules/03c_rseqc.snake"
include: "rules/03d_results_riboseq.snake"
# BW with deeptools bamCoverage
include: "rules/04a_bamCoverage_stranded.snake"
include: "rules/04b_bw_UCSC_URL_stranded.snake"
include: "rules/04c_riboseq_qc_plots.snake"
