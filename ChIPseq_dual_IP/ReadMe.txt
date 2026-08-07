# in project folder sync ChIPseq_dual_IP pipeline

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/ChIPseq_dual_IP/setup_pipeline.sh
# run script with one of the profiles: Bodhi | Alpine
bash setup_pipeline.sh Bodhi

# edit ChIPseq_dual_IP_samples.yaml with your sample information
# also review/edit UnStranded_matrix.yaml (computeMatrix/heatmap parameters, steps 5-7 below,
#   run as a second snakemake workflow against workflow/UnStranded_matrix.smk)
# run script by submitting to Bodhi (AMC) or Alpine (Boulder)

!!! Genome must point to ref yaml file with 2 genomes !!!
!!! see workflow/ref/hg38_polIII_SEQL_hg38.yaml for example !!!


### Pipeline Overview ###
1. Quality Control & Preprocessing

1a FASTQC — Quality control of raw data
1b bbtools clumpify — Remove optical/PCR duplicate sequences
1c bbtools bbduk — Remove Illumina adapters

2. Alignment & BAM Processing

2a bowtie2 — Align samples to two reference genomes in sequence; UNALIGNED=True routes only reads unmapped
   after the first pass into the second alignment (sequential/exclusive dual IP), UNALIGNED=False independently
   aligns all reads to both genomes
2b samtools — BAM filtering and sorting - both alignments
2c samtools (optional) — Subsample BAM files to normalize read depth across paired sample groups
   (NORM: "subsample") - both alignments
2d bam URLs (Bodhi option) — Copy BAM files to sandbox with URLs - both alignments

3. Feature Counting & Quality Metrics

3a subReads featureCounts — Count reads over annotated genes (protein-coding/scalefactor breakdown applied later in the summary report) - both alignments
3b deepTools fragmentSize — Report read and insert size distributions - both alignments
3c Rscript — Collect logs into summary report table - both alignments
3d MultiQC (optional) — Generate interim HTML summary report - both alignments

4. Coverage Track Generation - both alignments

4a deepTools bamCoverage — Create bigWig coverage files - both alignments (optional scalefactor norm of one alignment to the other)
4b bigWig URLs (Bodhi option) — Copy bigWig to sandbox and generate UCSC URLs - both alignments

5. Matrix Computation - primary genome

5a deepTools computeMatrix — Create matrix files - main alignment
5b BW URLs (Bodhi option) — Copy matrix to sandbox with URLs for bentools - main alignment

6. Visualization - primary genome

6a deepTools plotHeatmap — Generate plots and heatmaps from matrix files

7. Final Report

7a MultiQC — Create comprehensive HTML summary report
