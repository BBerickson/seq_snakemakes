#!/usr/bin/env bash

#BSUB -J RNAseq
#BSUB -o logs/RNAseq_%J.out
#BSUB -e logs/RNAseq_%J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4] span[hosts=1]"

set -o nounset -o pipefail -o errexit -x

# Create necessary directories
mkdir -p logs

# Configuration
PROFILE="workflow/profiles/Bodhi"
LSF_CONFIG="workflow/profiles/Bodhi/Bodhi_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="RNAseq_samples.yaml"
MATRIX_FILE="Stranded_matrix.yaml"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/star.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/Stranded_matrix.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} ${MATRIX_FILE} \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

