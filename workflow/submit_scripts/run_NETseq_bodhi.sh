#!/usr/bin/env bash

#BSUB -J NETseq
#BSUB -o logs/NETseq_%J.out
#BSUB -e logs/NETseq_%J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4] span[hosts=1]"

set -o nounset -o pipefail -o errexit -x

# Create necessary directories
mkdir -p logs

# Configuration
PROFILE="workflow/profiles/Bodhi"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="NETseq_samples.yaml"
MATRIX_FILE="Stranded_matrix.yaml"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/NETseq.smk \
    --configfile ${SAMPLES_FILE} \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/Stranded_matrix.smk \
    --configfile ${SAMPLES_FILE} ${MATRIX_FILE} \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

