# in project folder sync RNAseq pipeline

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/RNAseq/setup_pipeline.sh
# run script with one of the profiles: Bodhi | Alpine
bash setup_pipeline.sh Bodhi

# edit RNAseq_samples.yaml with your sample information
# run script by submitting to Bodhi (lsf) or Alpine (slerm) 

!!!STAR is not set up for mixed genome alignment!!!


### Pipeline Overview ###
1. Quality Control & Preprocessing

1a FASTQC — Quality control of raw data
1b bbtools clumpify — Remove optical/PCR duplicate sequences
1c bbtools bbduk — Remove Illumina adapters

2. Alignment & BAM Processing

2a HISAT2/STAR — Align samples to reference genome
2b samtools — BAM filtering and sorting
2c samtools (optional) — Subsample BAM for small RNA masking
2d bam URLs (Bodhi option) — Copy BAM files to sandbox with URLs

3. Feature Counting & Quality Metrics

3a subReads featureCounts — Count protein-coding genes
3b deepTools fragmentSize — Report read and insert size distributions
3c Rscript — Collect logs into summary report table
3d Rmarkdown (optional) — Generate HTML summary report

4. Coverage Track Generation

4a deepTools bamCoverage — Create stranded bigWig coverage files
4b bigWig URLs (Bodhi option) — Copy bigWig to sandbox and generate UCSC URLs

5. Matrix Computation

5a deepTools computeMatrix — Create sense/antisense stranded matrix files
5b BW URLs (Bodhi option) — Copy matrix to sandbox with URLs for bentools

6. Visualization

6a deepTools plotHeatmap — Generate plots and heatmaps from matrix files

7. Final Report

7a Rmarkdown — Create comprehensive HTML summary report
