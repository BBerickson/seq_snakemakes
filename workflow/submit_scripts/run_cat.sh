#!/usr/bin/env bash

#SBATCH --job-name=samples_cat
#SBATCH --output=logs/samples_cat_%j.out
#SBATCH --error=logs/samples_cat_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=normal 
#SBATCH --qos=normal

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

# Configuration
DATASET="samples_cat"  # Set your dataset name here

PROFILE="workflow/profiles/Bodhi_SLURM"
LSF_CONFIG="workflow/profiles/Bodhi_SLURM/Bodhi_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="${DATASET}.yaml"
SNAKE_FILE="${DATASET}.smk"
SINGULARITY_PREFIX="/projects/${USER}/.snakemake/singularity"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/${SNAKE_FILE} \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"



