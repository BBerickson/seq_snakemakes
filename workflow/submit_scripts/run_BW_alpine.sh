#!/usr/bin/env bash

#SBATCH --job-name=BW
#SBATCH --output=logs/BW_%j.out
#SBATCH --error=logs/BW_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

module purge

# Load modules
module load python/3.10.2
module load singularity/3.7.4
module load miniforge

conda activate bioinfo_env # snakemake

# Configuration
DATASET="ChIPseq"  # Set your dataset name here
MATRIXSET="UnStranded" # set the type here

PROFILE="workflow/profiles/Alpine"
LSF_CONFIG="workflow/profiles/Alpine/Alpine_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="${DATASET}_samples.yaml"
MATRIX_FILE="${MATRIXSET}_matrix.yaml"
MATRIX_SNAKE="${MATRIXSET}_matrix.smk"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/${MATRIX_SNAKE} ${LSF_CONFIG} \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} ${MATRIX_FILE} \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

