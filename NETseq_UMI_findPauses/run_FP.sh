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
module load R/4.2.2
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
        -R "rusage[mem={params.memory}] span[hosts=1]"
        -n {threads}
        -q rna '

    /beevol/home/erickson/.local/bin/snakemake \
        --snakefile $snake_file \
        --drmaa "$args" \
        --jobs 100 \
        --configfile $config_file
}

# Run pipeline to process UMI Netseq 
pipe_dir=pipelines
pipe_dir1=pipelines/ref

samples=samples.yaml
samples2=samples_FP.yaml

genome=$pipe_dir1/hg38_FP.yaml
snake=$pipe_dir/NETseq_pauses.snake
run_snakemake $snake "$samples $samples2 $genome"


