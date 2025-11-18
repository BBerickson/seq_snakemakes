#!/usr/bin/env bash

#BSUB -J samples_cat
#BSUB -o logs/samples_cat_%J.out
#BSUB -e logs/samples_cat_%J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4] span[hosts=1]"

set -o nounset -o pipefail -o errexit -x

# Shared Singularity cache for all Snakemake projects
SINGULARITY_PREFIX="/beevol/home/${USER}/.singularity_cache"

# Create necessary directories
mkdir -p logs

# Configuration
DATASET="samples_cat"  # Set your dataset name here

PROFILE="workflow/profiles/Bodhi"
LSF_CONFIG="workflow/profiles/Bodhi/Bodhi_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="${DATASET}.yaml"
SNAKE_FILE="${DATASET}.smk"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/${SNAKE_FILE} \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"



