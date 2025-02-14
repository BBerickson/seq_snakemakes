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

# Run pipeline to process ChIPseq reads
pipe_dir=pipelines
pipe_dir1=pipelines/ref
# index and configs
snake=$pipe_dir/NETseq.snake
genome=$pipe_dir1/hg38.yaml
samples=samples.yaml

run_snakemake $snake "$samples $genome"

# snake=$pipe_dir/Stranded_matrix_offset_nogroup.snake
# config=$pipe_dir/Stranded_matrix_UMI.yaml
# run_snakemake $snake "$samples $config $genome"

