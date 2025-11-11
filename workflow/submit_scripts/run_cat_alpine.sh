#!/usr/bin/env bash

#SBATCH --job-name=samples_cat
#SBATCH --output=logs/samples_cat_%j.out
#SBATCH --error=logs/samples_cat_%j.err
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
DATASET="samples_cat"  # Set your dataset name here

PROFILE="workflow/profiles/Alpine"
LSF_CONFIG="workflow/profiles/Alpine/Alpine_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="${DATASET}.yaml"
SNAKE_FILE="${DATASET}.smk"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/${SNAKE_FILE} \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"



