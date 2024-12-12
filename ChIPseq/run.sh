#!/usr/bin/env bash

#BSUB -J snake
#BSUB -o logs/snake_%J.out
#BSUB -e logs/snake_%J.err
#BSUB -q rna

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
        -n {threads}
        -q rna '

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
snake=$pipe_dir/ChIPseq_PE_spikeIN.snake
#snake=$pipe_dir/ChIPseq_PE.snake
#snake=$pipe_dir/ChIPseq_SE.snake
genome=$pipe_dir1/hg38_mm10.yaml
#genome=$pipe_dir1/hg38_scer3.yaml
#genome=$pipe_dir1/hg38.yaml
samples=samples.yaml

run_snakemake $snake "$samples $genome"

# Run pipeline to make table files of sample
snake=$pipe_dir/UnStranded_matrix.snake
config=$pipe_dir/UnStranded_matrix.yaml

run_snakemake $snake "$samples $config $genome"

# Run pipeline to make table files of spikeIN
# genome=$pipe_dir1/mm10.yaml
# snake=$pipe_dir/ChIPseq_matrix2.snake
# config=$pipe_dir/ChIPseq_matrix.yaml
# 
# run_snakemake $snake "$samples $config $genome"

