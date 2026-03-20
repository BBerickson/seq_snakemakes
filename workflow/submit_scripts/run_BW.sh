#!/usr/bin/env bash

#SBATCH --job-name=BW
#SBATCH --output=logs/BW_%j.out
#SBATCH --error=logs/BW_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=normal 
#SBATCH --qos=normal

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

# Configuration
DATASET="ChIPseq"  # Set your dataset name here
MATRIXSET="UnStranded" # set the type here

PROFILE="workflow/profiles/Bodhi_SLURM"
LSF_CONFIG="workflow/profiles/Bodhi_SLURM/Bodhi_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="${DATASET}_samples.yaml"
MATRIX_FILE="${MATRIXSET}_matrix.yaml"
MATRIX_SNAKE="${MATRIXSET}_matrix.smk"
SINGULARITY_PREFIX="/projects/${USER}/.snakemake/singularity"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/${MATRIX_SNAKE} ${LSF_CONFIG} \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} ${MATRIX_FILE} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

