# Process DNAseq and RNAseq, from raw data to bam to bigwig and deeptools matrix files

This snakemake pipeline runs in a container and is set up to run on AMC Bodhi (LSF) and Bolder Alpine (SLURM)

## Features

1. Quality Control & Preprocessing
  a FastQC on R1 and R2 (if available)
  b Clean up fastq files:
    -UMI workflow: extract UMI from R2/R1, use cutadapt to filter for adaptor-containing reads
    -OR Standard workflow:
      -clumpify: remove duplicates based on exact sequence/cluster location
      -bbduk: trim adaptors
  
2. align paired or single end fastqs in R1R2 or R2R1 orentation, mixed genomes and dual IP saported, outputs bam filtered file(s)
  a bowtie | hisat | star
    - splits out mixed genome alignment
    - sequenctual or independent dual IP alignment
  b filter bam files
    - samtools filtering
    - multimappers max k filtering
    - UMI filtering 
    - subset bams filtering
    - mask/blacklist regionsfiltering 
  c optional cp bams to amc-sandbox and makes a URL ref file
  
3. Statistics & QC reporting
  a Collect metrics: 
    - featureCounts
    - fragment size
    - preprocessing stats
    - alignment rates 
    - filtering stats
  b optional compile into summary table and generate QC R Markdown HTML report 

4. genome covrage bigwig files with deeptools, 
  - normalized: none, CPM, RPKM, scalefactor (dual IP), spikeIN +/- inputs (mixed genomes), chrM, featureCounts Geneid: ("snRNA","snoRNA","rRNA" ...)
  - unstranded | stranded pos neg files
  - optional cp to amc-sandbox and make URLs for UCSC browser
  
5. deeptools matrix files 
  a scale-regions: 543, reference-point: 5, 3, PI
  b optional cp to amc-sandbox and make URLs for bentools
  
6. Overview rmarkdown HTML QC report + plots
  a QC report 
  b deeptools profile, heatmap, cluster plots
    - merging and plotting grouped samples
    - pdf of plots per region
 
## Usage Guide

  Open directory of experement type and follow instructions in the ReadMe to copy from github the setup script to the project directory on bodhi or alpine
  Run the setup script that will copy all needed files for experement type to the project directory
  Edit sample.yaml and if needed the matrix.yaml
  Run the snakemake script
  
