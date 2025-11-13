# ===== Snake file for processing Bowtie ================================

# Configure shell for all rules
shell.executable("/bin/bash")
shell.prefix("source ~/.bash_profile; set -o nounset -o pipefail -o errexit -x; ")

# python packages
import os
import sys
import yaml
from pathlib import Path

# Include custom Python functions
include: workflow.source_path("scripts/funs.py")

# ------------------------------------------------------------------------------
# Assign parameters from configs
# ------------------------------------------------------------------------------

# From main config
PROJ         = config.get("PROJ")
RAW_DATA     = config.get("RAW_DATA")
ALL_SAMPLES  = config.get("SAMPLES")
SEQ_DATE     = config.get("SEQ_DATE")
PAIRED       = config.get("PAIRED")

# Directories for data and scripts
FASTQ_DIR = PROJ + "/raw_fastq_links"

os.makedirs(FASTQ_DIR, exist_ok = True)

def process_samples(all_samples):
    SAMPLES = {}    # {sample_name: [fastq_files]}
    
    # Process the nested structure: SAMPLES -> SECTION -> PAIR -> SAMPLE -> fastq_files
    for section_name, section_data in all_samples.items():
        for pair_name, pair_data in section_data.items():
            for sample_name, fastq_files in pair_data.items():
                # Each sample_name maps to a list of fastq files
                SAMPLES[sample_name] = fastq_files
    
    return SAMPLES

def create_fastq_groups(sample_groups, raw_data, fastq_dir, paired):
    """
    Transform sample_groups into a dictionary with R1 and R2 file lists
    Returns: {group: {"R1": [list_of_R1_files], "R2": [list_of_R2_files]}}
    """
    fastq_groups = {}
    
    for group, samples in sample_groups.items():
        r1_files = []
        r2_files = []
        
        for sample in samples:
            fastqs = _get_fqs(sample, raw_data, fastq_dir, paired)
            r1_files.append(fastqs[0])  # Always has R1
            if len(fastqs) > 1:  # Has R2
                r2_files.append(fastqs[1])
        
        fastq_groups[group] = {
            "R1": r1_files,
            "R2": r2_files if r2_files else None  # None if no R2 files
        }
    
    return fastq_groups

# Simplify ALL_SAMPLES dictionary
SAMPLES = process_samples(ALL_SAMPLES)

# Create the processed dictionary
FASTQ_GROUPS = create_fastq_groups(SAMPLES, RAW_DATA, FASTQ_DIR, PAIRED)

# Print summary of samples and groups
print("SAMPLES (%s): %s\n" % (len(SAMPLES), SAMPLES))
print("FASTQ_GROUPS (%s): %s\n" % (len(FASTQ_GROUPS), FASTQ_GROUPS))

# Wildcard constraints
WILDCARD_REGEX = r"[a-zA-Z0-9_\-]+" # Matches alphanumeric characters, underscores, and hyphens

wildcard_constraints:
    sample = WILDCARD_REGEX,
    newnam = WILDCARD_REGEX


# Final output files
def get_expected_outputs():
    outputs = []
    for newnam, files in FASTQ_GROUPS.items():
        outputs.append(f"fastqs/{newnam}_R1_001.fastq.gz")
        if files["R2"]:  # Only add R2 if it exists
            outputs.append(f"fastqs/{newnam}_R2_001.fastq.gz")
    return outputs

rule all:
    input:
        get_expected_outputs()


# Run concatenate
include: "rules/00_concatenate.snake"

