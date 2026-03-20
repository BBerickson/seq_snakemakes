#!/usr/bin/env bash

#SBATCH --job-name=RNAseq
#SBATCH --output=logs/RNAseq_%j.out
#SBATCH --error=logs/RNAseq_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=normal 
#SBATCH --qos=normal

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

# Configuration
PROFILE="workflow/profiles/Bodhi_SLURM"
LSF_CONFIG="workflow/profiles/Bodhi_SLURM/Bodhi_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="RNAseq_samples.yaml"
MATRIX_FILE="Stranded_matrix.yaml"
SINGULARITY_PREFIX="/projects/${USER}/.snakemake/singularity"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/RNAseq.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/Stranded_matrix.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG}  ${MATRIX_FILE} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

