# Run on LSF system with snakemake in a singularity container

#1 Quality Control & Preprocessing
  # FastQC on R1 and R2 (if available)
  # Clean up fastq files:
    # UMI workflow: extract UMI from R2/R1, use cutadapt to filter for adaptor-containing reads
    # OR Standard workflow:
      # clumpify: remove duplicates based on exact sequence/cluster location
      # bbduk: trim adaptors
  
#2 align paired or single end fastqs, R1R2 or R2R1, clean up and outputs bam file(s) (mask is only sample not spikeIN)
  #a bowtie | hisat | star
  #b clean up, filter, sort bam files
    # +/- UMI filtering, +/- subset bams, +/- mask blacklist regions, filter/split genomes -> output bam(s)
  #c optional cp bams to amc-sandbox and makes a URL ref file
  
#3 Statistics & QC reporting
  # Collect metrics: featureCounts, fragment size, PCA, preprocessing stats, alignment rates, filtering stats
  # Compile into summary table
  # Generate QC R Markdown HTML report 

#4 genome covrage bigwig files, 
  # normalized: none, CPM, RPKM, spikeIN +/- inputs, chrM, 
  # unstranded | stranded pos neg files
  # cp to amc-sandbox and make URLs for UCSC browser
  
#5 deeptools matrix files 
  # scale-regions: 543, reference-point: 5, 3, PI
  # cp to amc-sandbox and make URLs for bentools
  
#6 deeptools profile, heatmap, cluster plots
  # merging and plotting grouped samples
  # pdf of plots per region
#7 rmarkdown HTML QC report + plots