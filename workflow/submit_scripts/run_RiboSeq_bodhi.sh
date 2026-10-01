#!/usr/bin/env bash

#SBATCH --job-name=RiboSeq
#SBATCH --output=logs/RiboSeq_%j.out
#SBATCH --error=logs/RiboSeq_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=normal
#SBATCH --qos=normal

set -o nounset -o pipefail -o errexit -x

# Shared Singularity cache for all Snakemake projects
SINGULARITY_PREFIX="/beevol/home/${USER}/.singularity_cache"

# Create necessary directories
mkdir -p logs

# Configuration
PROFILE="workflow/profiles/Bodhi"
LSF_CONFIG="workflow/profiles/Bodhi/Bodhi_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="RiboSeq_samples.yaml"

# CONTAINER is normally set at the profile level (Bodhi_config.yaml), pointing
# at the base bioinfopipe image -- Ribo-seq needs its own dedicated image
# (Trim Galore/SortMeRNA/Salmon/RSeQC/riboWaltz/anota2seq/xtail), so override
# it here rather than duplicating an entire profile directory just to change
# one image reference. Update this to wherever you push the Ribo-seq image.
RIBOSEQ_CONTAINER="docker://bberickson/riboseq:1.0.0"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/RiboSeq.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}" CONTAINER="${RIBOSEQ_CONTAINER}"
