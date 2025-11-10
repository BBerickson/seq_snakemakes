#!/usr/bin/env bash

#BSUB -J ChIPseq
#BSUB -o logs/ChIPseq_%J.out
#BSUB -e logs/ChIPseq_%J.err
#BSUB -n 1
#BSUB -R "rusage[mem=4] span[hosts=1]"

set -o nounset -o pipefail -o errexit -x

# Create necessary directories
mkdir -p logs

# Configuration
PROFILE="workflow/profiles/Bodhi"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="samples_cat.yaml"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/samples_cat.smk \
    --configfile ${SAMPLES_FILE} \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"



