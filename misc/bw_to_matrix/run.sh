#!/usr/bin/env bash

#BSUB -J snake
#BSUB -o logs/snake_%J.out
#BSUB -e logs/snake_%J.err

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

# Run pipeline to process bigwig files to matrix files
# index and configs
samples=samples.yaml

# Run pipeline to make table files of sample
# snake=pipelines/BW_matrix.snake
# config=pipelines/UnStranded_matrix.yaml
# snake=pipelines/BW_bidirectional_matrix.snake
# config=pipelines/bidirectional_matrix.yaml
snake=pipelines/BW_Stranded_matrix.snake
config=pipelines/Stranded_matrix.yaml

run_snakemake $snake "$samples $config"


