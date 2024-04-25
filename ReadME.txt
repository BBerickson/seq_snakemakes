# pipeline notes

#1 fastQC R1 and R2 if avalible

#1 clean up fastq
  #UMI: extract UMI from read 2, using cutadapt look for and keep reads with adaptors
  #a clumpify: remove duplicates based on exact sequence and location of cluster, other options
  #b bbduk: trim adaptors, provided or looks at fastq header and gets adaptor from list
  
#2 align and fastqs, cleans up and outputs bam final form(s) (mask is only sample not spikeIN)
  #a bowtie, hisat, star
  #b clean up, filter, sort bam
    # filter/split based on common chromosomes -> output bam(s)
    # mask filter/split based on common chromosomes -> output bam(s) 
    # UMI filtering, +/- subset bams, filter/split -> output bam(s)
    # UMI filtering, +/- subset bams, mask filter/split -> output bam(s)
  #c cp bams to sandbox and makes a URL
  
#3 pull together all log files into a results file
  # used to indicate if there is a spike in file to separate out or not
  # UMI processing produces diffrnet logs so its own rule is needed

#4 make bigwig files from bams, normalizing, cp to sandbox and make URLs for UCSC browser
  # +/- stranded 
  # UMI can't be renamed with group name so nogroup rule
  
#5 make matrix files from bigwigs, cp to sandbox and make URLs for bentools
  # scale-regions or reference-point