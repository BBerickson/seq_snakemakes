#!/usr/bin/env bash

#SBATCH --job-name=RiboSeq_TE
#SBATCH --output=logs/RiboSeq_TE_%j.out
#SBATCH --error=logs/RiboSeq_TE_%j.err
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=normal
#SBATCH --qos=normal

set -o nounset -o pipefail -o errexit -x

# Shared Singularity cache for all Snakemake projects
SINGULARITY_PREFIX="/beevol/home/${USER}/.singularity_cache"

mkdir -p logs

# Configuration -- run RiboSeq.smk/RiboSeq_UMI.smk and RNAseq.smk to
# completion FIRST; this pipeline only joins their already-completed outputs.
PROFILE="workflow/profiles/Bodhi"
LSF_CONFIG="workflow/profiles/Bodhi/Bodhi_config.yaml"
SAMPLES_FILE="RiboSeq_TE_samples.yaml"

# Same dedicated Ribo-seq image as run_RiboSeq_bodhi.sh (needs anota2seq/xtail,
# and clusterProfiler + org.*.eg.db if run_go_enrichment is turned on).
RIBOSEQ_CONTAINER="docker://bberickson/riboseq:1.0.0"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/RiboSeq_TE.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config CONTAINER="${RIBOSEQ_CONTAINER}"
