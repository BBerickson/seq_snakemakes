#!/usr/bin/env bash

#BSUB -J ChIPseq
#BSUB -o logs/ChIPseq_%J.out
#BSUB -e logs/ChIPseq_%J.err

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


# index and configs
snake=pipelines/ChIPseq_duel_IP.snake
samples=samples.yaml

# Run pipeline to process ChIPseq reads
run_snakemake $snake "$samples"

# configs 
snake=pipelines/UnStranded_matrix.snake
config=pipelines/UnStranded_matrix.yaml

# Run pipeline to make table files of sample
run_snakemake $snake "$samples $config"
