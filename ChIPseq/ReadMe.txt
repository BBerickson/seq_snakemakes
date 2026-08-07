# in project folder get ChIPseq setup script

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/ChIPseq/setup_pipeline.sh
# run script with one of the profiles: Bodhi | Alpine
bash setup_pipeline.sh Bodhi

# edit ChIPseq_samples.yaml with your sample information
# also review/edit UnStranded_matrix.yaml (computeMatrix/heatmap parameters, steps 5-7 below,
#   run as a second snakemake workflow against workflow/UnStranded_matrix.smk)
# run script by submitting to Bodhi (AMC) or Alpine (Boulder)

!!! Mixed genomes must point to ref yaml file with 2 genomes that points to 2 ref yamls !!!
!!! see workflow/ref/hg38_mm39.yaml for example !!!

### Pipeline Overview ###
1. Quality Control & Preprocessing

1a FASTQC — Quality control of raw data
1b bbtools clumpify — Remove optical/PCR duplicate sequences
1c bbtools bbduk — Remove Illumina adapters

2. Alignment & BAM Processing

2a bowtie2 — Align samples to reference genome, or mixed reference genomes 
2b samtools — BAM filtering and sorting, (spliting reference genomes)
2c samtools (optional) — Subsample aligned BAM files to normalize read depth across paired sample groups (NORM: "subsample")
2d bam URLs (Bodhi option) — Copy BAM files to sandbox with URLs

3. Feature Counting & Quality Metrics

3a subReads featureCounts — Count reads over annotated genes (protein-coding/scalefactor breakdown applied later in the summary report)
3b deepTools fragmentSize — Report read and insert size distributions
3c Rscript — Collect logs into summary report table
3d MultiQC (optional) — Generate interim HTML summary report

4. Coverage Track Generation - all configured genomes (e.g. spike-in indexes get auto CPM-normalized tracks)

4a deepTools bamCoverage — Create bigWig coverage files 
4b bigWig URLs (Bodhi option) — Copy bigWig to sandbox and generate UCSC URLs

5. Matrix Computation - primary genome

5a deepTools computeMatrix — Create matrix files
5b BW URLs (Bodhi option) — Copy matrix to sandbox with URLs for bentools

6. Visualization - primary genome

6a deepTools plotHeatmap — Generate plots and heatmaps from matrix files

7. Final Report

7a MultiQC — Create comprehensive HTML summary report
