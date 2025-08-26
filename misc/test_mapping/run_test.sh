#!/usr/bin/env bash

#BSUB -J snake
#BSUB -o logs/snake_test_%J.out
#BSUB -e logs/snake_test_%J.err


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

# Run pipeline to process reads
# index and configs
snake=pipelines/bowtie_test.snake
#snake=pipelines/NETseq_UMI_test.snake
samples=samples.yaml

run_snakemake $snake "$samples"


