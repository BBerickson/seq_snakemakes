#!/usr/bin/env bash

#SBATCH --job-name=rMATS
#SBATCH --output=logs/rmats_%j.out
#SBATCH --error=logs/rmats_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=normal 
#SBATCH --qos=normal

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

# Configuration
DATASET="RNAseq"  # Set your dataset name here

PROFILE="workflow/profiles/Bodhi_SLURM"
LSF_CONFIG="workflow/profiles/Bodhi_SLURM/Bodhi_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="${DATASET}_samples.yaml"
SAMPLES_RMATS="samples_rMATS.yaml"
SINGULARITY_PREFIX="/projects/${USER}/.snakemake/singularity"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/rmats.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} ${SAMPLES_RMATS} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

