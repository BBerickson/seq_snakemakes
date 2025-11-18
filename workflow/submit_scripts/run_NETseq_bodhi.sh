#!/usr/bin/env bash

#BSUB -J NETseq
#BSUB -o logs/NETseq_%J.out
#BSUB -e logs/NETseq_%J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4] span[hosts=1]"

set -o nounset -o pipefail -o errexit -x

# Shared Singularity cache for all Snakemake projects
SINGULARITY_PREFIX="/beevol/home/${USER}/.singularity_cache"

# Create necessary directories
mkdir -p logs
mkdir -p "${SINGULARITY_CACHEDIR}" 

# Configuration
PROFILE="workflow/profiles/Bodhi"
LSF_CONFIG="workflow/profiles/Bodhi/Bodhi_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="NETseq_samples.yaml"
MATRIX_FILE="Stranded_matrix.yaml"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/NETseq.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/Stranded_matrix.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} ${MATRIX_FILE} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

