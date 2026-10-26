# in project folder sync RiboSeq_UMI pipeline

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/RiboSeq_UMI/setup_pipeline.sh
# run script with one of the profiles: Bodhi
bash setup_pipeline.sh Bodhi

# edit RiboSeq_UMI_samples.yaml with your sample information (UMI pattern,
# extraction orientation) and run run_RiboSeq_UMI_bodhi.sh by submitting to
# Bodhi (AMC)

# For a non-UMI prep (classic Ingolia-style circularization) instead, use
# the separate RiboSeq pipeline (its own setup_pipeline.sh, same repo).


### Pipeline Overview -- Phase 1: Ribo-seq alignment + QC (UMI variant) ###

1. Quality Control & Preprocessing

1a FASTQC — Quality control of raw data
1a' umi_tools extract — Extracts UMI into the read name before trimming
    (UMI-specific step; same umi_pattern/ORIENTATION convention as the
    other UMI pipelines in this repo, e.g. BRUseq_UMI)
1b cutadapt — Adapter/quality trimming (on the UMI-extracted fastq)
1c bbsplit (optional) — Contaminant removal (phiX etc.), skip_bbsplit: true to disable
1d SortMeRNA — rRNA removal (chunked, rebuilds its index once per chunk)
1e Salmon (advisory only) — Infers library strandedness; reports to
    PROJ/stats/PROJ_strandedness.tsv, doesn't self-configure ORIENTATION

2. Alignment & BAM Processing

2a STAR (2-pass) — Aligns to genome AND transcriptome
   (--quantMode TranscriptomeSAM, for riboWaltz below)
2b umi_tools dedup — PCR/optical duplicate removal by UMI (genome BAM;
   UMI-specific step, no equivalent in the non-UMI RiboSeq pipeline)
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
3e Rscript — Collect logs into summary report table (UMI dedup stats
   included, unlike the non-UMI variant)
3f P-site-restricted, CDS-only gene counts (optional, enable_psite_cds_counts:
   true) — riboWaltz's cds_coverage(): counts only in-frame P-sites inside
   the annotated CDS, a stricter alternative to 3a's whole-exon counts for
   downstream TE analysis. Output: PROJ/counts/{sample}_..._psite_cds_
   featureCounts.tsv, same column schema as 3a's output (drop-in for
   Phase 2's CMD_PARAMS.use_psite_counts toggle). Expect substantially
   lower totals than 3a -- that's CDS/in-frame restriction + multi-mapping
   resolution working as intended, not a bug; cross-check against 3d's
   intron % and 3c's in_frame_pct if the magnitude looks surprising.
3g Differential ribosome occupancy (optional, CONTRASTS set in
   RiboSeq_UMI_samples.yaml) — plain two-group DESeq2 test on this
   project's OWN featureCounts, by SAMPLES nickname, no paired RNA-seq
   needed or used. NOT translational efficiency -- without RNA-seq there's
   no way to tell "more ribosomes on this transcript" apart from "more of
   this transcript exists" (see Phase 2 below for the real thing) -- but a
   real, standard thing to report when no paired RNA-seq run exists, with
   that caveat. Uses 3a's counts by default, or 3f's P-site/CDS-restricted
   counts via CMD_PARAMS.differential_occupancy.use_psite_counts.
3h MultiQC (optional) — Generate interim HTML summary report, including
   3g's MA/volcano plots if CONTRASTS is set

4. Coverage Track Generation

4a deepTools bamCoverage — Create stranded bigWig coverage files
   (single-end footprints, still stranded -- reflects transcript strand)
4b bigWig URLs (Bodhi option) — Copy bigWig to sandbox and generate UCSC URLs
4c deepTools multiBamSummary/plotPCA/plotCorrelation — Sample QC
   (replicate clustering/correlation)

5. Final Report

5a MultiQC (optional) — Create comprehensive HTML summary report


### Pipeline Overview -- Phase 2: Translational Efficiency (optional) ###

Joins a COMPLETED Ribo-seq run with a COMPLETED, PAIRED RNA-seq run into
anota2seq (primary, needs >=3 replicates/group) + xtail (secondary, no
replicate minimum) translational-efficiency analysis. Does no
alignment/quantification of its own -- only reads both runs' completed
featureCounts outputs. Identical to the non-UMI RiboSeq pipeline's Phase 2
-- RiboSeq_TE.smk doesn't care whether the Ribo-seq side used UMIs.

If you already know up front you want this, use setup_TE_pipeline.sh
INSTEAD of setup_pipeline.sh above, in a fresh/empty project folder -- it
sets up all three arms together in one call:

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/RiboSeq_UMI/setup_TE_pipeline.sh
bash setup_TE_pipeline.sh Bodhi

This creates riboseq/ and rnaseq/ subfolders (each a complete, self-
contained Ribo-seq/RNA-seq setup), with the TE join (RiboSeq_TE.smk) at
this top level. Separate subfolders for riboseq/rnaseq, not this same
folder -- Snakemake locks its whole working directory for a run, so two
snakemake invocations sharing one CWD would collide the moment both run
at once, which is exactly what the pipeline script below does. Then:
  1. cd riboseq, edit RiboSeq_UMI_samples.yaml with your sample information
     (UMI pattern, extraction orientation), and run run_RiboSeq_UMI_bodhi.sh.
  2. cd rnaseq, edit RNAseq_samples.yaml with your paired RNA-seq sample
     information, and run one of its submit scripts (run_RNAseq_Hisat_bodhi.sh /
     run_RNAseq_STAR_bodhi.sh).
  3. Once both arms are done, edit RiboSeq_TE_samples.yaml here (RIBO_PROJ
     -- riboseq/<its PROJ> -- and RNA_PROJ -- rnaseq/<its PROJ> --,
     CONTRASTS by sample nickname, CMD_PARAMS.use_psite_counts if using
     step 3f's counts instead of 3a's), then either:
       - run run_RiboSeq_TE_bodhi.sh, once both arms are already done, or
       - run ./run_RiboSeq_TE_pipeline_bodhi.sh --ribo ./riboseq
         --rna ./rnaseq --te . --ribo-script run_RiboSeq_UMI_bodhi.sh to
         launch all three stages in one command -- the two independent
         arms run in parallel in their own subfolders, TE is deferred via
         a SLURM job dependency until both finish successfully.
         --ribo-script must be set here -- the orchestrator script
         defaults to the non-UMI run_RiboSeq_bodhi.sh. Note: --rna/--te
         are each required for that arm to be submitted at all -- passing
         only --ribo silently submits just the Ribo-seq arm, which is not
         a bug.

If you already have a completed Ribo-seq run from this folder's plain
setup_pipeline.sh and only decided you want TE afterward, don't run
setup_TE_pipeline.sh here -- it would set up a second, unrelated Ribo-seq
copy in riboseq/. Instead, get just RiboSeq_TE.smk + RNAseq.smk (e.g. from
a fresh setup_TE_pipeline.sh run elsewhere) and point RiboSeq_TE_samples.yaml's
RIBO_PROJ directly at this folder's existing run -- any path works.
