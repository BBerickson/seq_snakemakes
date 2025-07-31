#!/usr/bin/env bash

#BSUB -J snake
#BSUB -o logs/snake_cat_%J.out
#BSUB -e logs/snake_cat_%J.err


set -o nounset -o pipefail -o errexit -x

mkdir -p logs

# Load modules
. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load python/3.8.5

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
        --snakefile $snake_file \
        --drmaa "$args" \
        --jobs 100 \
        --configfile $config_file
}

# Run pipeline to process ChIPseq reads
pipe_dir=pipelines
# index and configs
snake=$pipe_dir/samples_cat.snake
samples=samples_cat.yaml

run_snakemake $snake "$samples"


