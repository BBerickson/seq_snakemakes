# in project folder get starting from Bam setup scripts


wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/misc/start_bam/setup_stranded_pipeline.sh
# or
wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/misc/start_bam/setup_Unstranded_pipeline.sh
# run script with one of the profiles: Bodhi | Alpine
bash setup_stranded_pipeline.sh Bodhi
or
bash setup_Unstranded_pipeline.sh Bodhi

# edit Bam_Stranded_samples.yaml (stranded) or Bam_UnStranded_samples.yaml (unstranded) with your sample information
# run script by submitting to Bodhi (AMC) or Alpine (Boulder)

### Pipeline Overview ###
1. obtaining bam and bai files

1. get_bams

2. BAM Processing

2b samtools — BAM filtering and sorting, (spliting reference genomes)
2c samtools (optional) — Subsample BAM for group-level read-depth normalization (stranded variant also masks abundant/filtered RNA loci first)
2d bam URLs (Bodhi option) — Copy BAM files to sandbox with URLs

3. Feature Counting & Quality Metrics

3a subReads featureCounts — Count reads over annotated genes (protein-coding/scalefactor breakdown applied later in the summary report)
3b deepTools fragmentSize — Report read and insert size distributions

4. Coverage Track Generation - primary genome

4a deepTools bamCoverage — Create bigWig coverage files 
4b bigWig URLs (Bodhi option) — Copy bigWig to sandbox and generate UCSC URLs

5. Matrix Computation - primary genome

5a deepTools computeMatrix — Create matrix files
5b BW URLs (Bodhi option) — Copy matrix to sandbox with URLs for bentools

6. Visualization - primary genome

6a deepTools plotHeatmap — Generate plots and heatmaps from matrix files

