#!/usr/bin/env bash

#BSUB -J snake
#BSUB -o logs/snake_%J.out
#BSUB -e logs/snake_%J.err


set -o nounset -o pipefail -o errexit -x

mkdir -p logs

# Load modules
. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load python/3.8.5
module load samtools/1.9
module load bbtools/39.01
module load bowtie2/2.3.2
module load R/4.3.3
module load fastqc/0.11.9
module load subread

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

# Run pipeline to process BrUseq reads
pipe_dir=pipelines
# index and configs
# snake=$pipe_dir/BRUseq_spikeIN_UMI.snake
snake=$pipe_dir/BRUseq_UMI.snake
samples=samples.yaml

run_snakemake $snake "$samples 

snake=$pipe_dir/Stranded_matrix.snake
config=$pipe_dir/Stranded_matrix.yaml
run_snakemake $snake "$samples $config 

# # Run pipeline to make table files of spikeIN

# snake=$pipe_dir/Stranded_matrix.snake
# config=$pipe_dir/Stranded_matrix_BRUseq.yaml
# 
# run_snakemake $snake "$samples $config"

