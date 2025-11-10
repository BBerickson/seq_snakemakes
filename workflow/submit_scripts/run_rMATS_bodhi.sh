#!/usr/bin/env bash

#BSUB -J rmats
#BSUB -o logs/rmats_%J.out
#BSUB -e logs/rmats_%J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4] span[hosts=1]"

set -o nounset -o pipefail -o errexit -x

# Create necessary directories
mkdir -p logs

# Configuration
DATASET="RNAseq"  # Set your dataset name here
MATRIXSET="UnStranded" # set the type here

PROFILE="workflow/profiles/Bodhi"
LSF_CONFIG="workflow/profiles/Bodhi/Bodhi_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="${DATASET}_samples.yaml"
SAMPLES_RMATS="samples_rMATS.yaml"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/rmats.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} ${SAMPLES_RMATS} \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

