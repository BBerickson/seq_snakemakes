#!/usr/bin/env bash

#SBATCH --job-name=ChIPseq
#SBATCH --output=logs/ChIPseq_%j.out
#SBATCH --error=logs/ChIPseq_%j.err
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
SAMPLES_FILE="ChIPseq_samples.yaml"
MATRIX_FILE="UnStranded_matrix.yaml"
SINGULARITY_PREFIX="/projects/${USER}/.snakemake/singularity"


# Run ChIPseq pipeline
snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/ChIPseq.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

# Run matrix generation pipeline
snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/UnStranded_matrix.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} ${MATRIX_FILE} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

