# in project folder sync RiboSeq pipeline

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/RiboSeq/setup_pipeline.sh
# run script with one of the profiles: Bodhi
bash setup_pipeline.sh Bodhi

# edit RiboSeq_samples.yaml with your sample information
# run run_RiboSeq_bodhi.sh by submitting to Bodhi (AMC)

# Ribo-seq is single-end, no UMIs (classic Ingolia-style circularization
# prep). For a UMI-based prep instead, use the separate RiboSeq_UMI
# pipeline (its own setup_pipeline.sh, same repo).


### Pipeline Overview -- Phase 1: Ribo-seq alignment + QC ###

1. Quality Control & Preprocessing

1a FASTQC — Quality control of raw data
1b cutadapt — Adapter/quality trimming
1c bbsplit (optional) — Contaminant removal (phiX etc.), skip_bbsplit: true to disable
1d SortMeRNA — rRNA removal (chunked, rebuilds its index once per chunk)
1e Salmon (advisory only) — Infers library strandedness; reports to
    PROJ/stats/PROJ_strandedness.tsv, doesn't self-configure ORIENTATION

2. Alignment & BAM Processing

2a STAR (2-pass) — Aligns to genome AND transcriptome
   (--quantMode TranscriptomeSAM, for riboWaltz below)
2b samtools — BAM filtering and sorting (genome BAM)
2c bam URLs (Bodhi option) — Copy BAM files to sandbox with URLs
2d multimappers2.py — Resolves multi-mapping on the TRANSCRIPTOME BAM
   specifically (keeps exactly 1 alignment/read) -- only runs if
   enable_psite_cds_counts: true (step 3f below); MAPQ does NOT reliably
   encode multi-mapping on this BAM, confirmed on real data, so this is a
   read-name-grouped filter, not a MAPQ threshold

3. Feature Counting & Quality Metrics

3a subReads featureCounts — Count reads over annotated genes (whole-gene
   exon overlap; protein-coding/scalefactor breakdown in the summary report)
3b deepTools fragmentSize — Footprint length distribution (~15-50nt)
3c riboWaltz — P-site offset calibration, periodicity, metagene profile,
   footprint length distribution (one summary row/sample: psite_offset,
   in_frame_pct, dominant_footprint_length)
3d RSeQC read_distribution.py — CDS/UTR/intron read-region breakdown
   (independent cross-check of library quality -- high intron %, if seen,
   reflects real rRNA/pre-mRNA contamination, not a pipeline issue)
3e Rscript — Collect logs into summary report table
3f P-site-restricted, CDS-only gene counts (optional, enable_psite_cds_counts:
   true) — riboWaltz's cds_coverage(): counts only in-frame P-sites inside
   the annotated CDS, a stricter alternative to 3a's whole-exon counts for
   downstream TE analysis. Output: PROJ/counts/{sample}_..._psite_cds_
   featureCounts.tsv, same column schema as 3a's output (drop-in for
   Phase 2's CMD_PARAMS.use_psite_counts toggle). Expect substantially
   lower totals than 3a -- that's CDS/in-frame restriction + multi-mapping
   resolution working as intended, not a bug; cross-check against 3d's
   intron % and 3c's in_frame_pct if the magnitude looks surprising.
3g MultiQC (optional) — Generate interim HTML summary report

4. Coverage Track Generation

4a deepTools bamCoverage — Create stranded bigWig coverage files
   (single-end footprints, still stranded -- reflects transcript strand)
4b bigWig URLs (Bodhi option) — Copy bigWig to sandbox and generate UCSC URLs
4c deepTools multiBamSummary/plotPCA/plotCorrelation — Sample QC
   (replicate clustering/correlation)

5. Final Report

5a MultiQC (optional) — Create comprehensive HTML summary report


### Pipeline Overview -- Phase 2: Translational Efficiency (optional) ###

Joins this COMPLETED Ribo-seq run with a COMPLETED, PAIRED RNA-seq run into
anota2seq (primary, needs >=3 replicates/group) + xtail (secondary, no
replicate minimum) translational-efficiency analysis. Does no
alignment/quantification of its own -- only reads both runs' completed
featureCounts outputs.

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/RiboSeq/setup_TE_pipeline.sh
bash setup_TE_pipeline.sh Bodhi

This creates sibling rnaseq/ and te/ subfolders (NOT this same folder --
Snakemake locks its whole working directory for a run, so the RNA-seq and
Ribo-seq snakemake invocations would collide the moment both run at once,
which is exactly what step 2's pipeline script below does), pulling the
paired RNA-seq arm (RNAseq.smk/star.smk) into rnaseq/ and the TE join
(RiboSeq_TE.smk) into te/. Then:
  1. cd rnaseq, edit RNAseq_samples.yaml with your paired RNA-seq sample
     information, and run one of its submit scripts (run_RNAseq_Hisat_bodhi.sh /
     run_RNAseq_STAR_bodhi.sh).
  2. Once both this Ribo-seq run and that RNA-seq run are done, cd te and
     edit RiboSeq_TE_samples.yaml (RIBO_PROJ -- one level up, this Ribo-seq
     run's own PROJ -- and RNA_PROJ -- ../rnaseq/<its PROJ> --, CONTRASTS by
     sample nickname, CMD_PARAMS.use_psite_counts if using step 3f's counts
     instead of 3a's), then either:
       - run run_RiboSeq_TE_bodhi.sh from inside te/, once both arms are
         already done, or
       - from this top folder, run
         ./run_RiboSeq_TE_pipeline_bodhi.sh --ribo . --rna ./rnaseq --te ./te
         to launch all three stages (this Ribo-seq run + the RNA-seq run +
         the TE join, each in its own subfolder) in one command -- the two
         independent arms run in parallel, TE is deferred via a SLURM job
         dependency until both finish successfully. Note: --rna/--te are
         each required for that arm to be submitted at all -- passing only
         --ribo (as you might by habit) silently submits just the Ribo-seq
         arm, which is not a bug.
