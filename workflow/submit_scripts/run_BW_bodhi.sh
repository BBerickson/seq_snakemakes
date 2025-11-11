#!/usr/bin/env bash

#BSUB -J BW
#BSUB -o logs/BW_%J.out
#BSUB -e logs/BW_%J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4] span[hosts=1]"

set -o nounset -o pipefail -o errexit -x

# Create necessary directories
mkdir -p logs

# Configuration
DATASET="ChIPseq"  # Set your dataset name here
MATRIXSET="UnStranded" # set the type here

PROFILE="workflow/profiles/Bodhi"
LSF_CONFIG="workflow/profiles/Bodhi/Bodhi_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="${DATASET}_samples.yaml"
MATRIX_FILE="${MATRIXSET}_matrix.yaml"
MATRIX_SNAKE="${MATRIXSET}_matrix.smk"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/${MATRIX_SNAKE} ${LSF_CONFIG} \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} ${MATRIX_FILE} \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

