#!/usr/bin/env bash

#BSUB -J BRUseq
#BSUB -o logs/BRUseq_%J.out
#BSUB -e logs/BRUseq_%J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4] span[hosts=1]"

set -o nounset -o pipefail -o errexit -x

# Create necessary directories
mkdir -p logs

# Configuration
PROFILE="profiles/Bodhi"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="BRUseq_samples.yaml"
MATRIX_FILE="Stranded_matrix.yaml"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/BRUseq.smk \
    --configfile ${SAMPLES_FILE} \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/Stranded_matrix.smk \
    --configfile ${SAMPLES_FILE} ${MATRIX_FILE} \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

