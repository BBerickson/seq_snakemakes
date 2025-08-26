#!/usr/bin/env bash

#BSUB -J rmats
#BSUB -o logs/rmats_%J.out
#BSUB -e logs/rmats_%J.err


set -o nounset -o pipefail -o errexit -x

mkdir -p logs

bind_dir='/beevol/home'
ssh_key_dir='$HOME/.ssh'

# Function to run snakemake
run_snakemake() {
    local snake_file=$1
    local config_file=$2

    args='
        -oo {log.out} 
        -eo {log.err} 
        -J {params.job_name}
        -R "rusage[mem={resources.memory}] span[hosts=1]"
        -n {threads} '

    snakemake \
        --use-singularity \
        --singularity-args "--bind $bind_dir" \
        --snakefile $snake_file \
        --drmaa "$args" \
        --jobs 100 \
        --config SSH_KEY_DIR="$ssh_key_dir" \
        --configfile $config_file
}

# Run pipeline to process RNAseq reads
# index and configs
snake=pipelines/rmats.snake
samples=samples.yaml
samples2=samples_rMATS.yaml

run_snakemake $snake "$samples $samples2"


