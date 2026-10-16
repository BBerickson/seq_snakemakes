#!/usr/bin/env Rscript
# P-site-restricted, CDS-only gene-level counting -- a stricter alternative
# to plain featureCounts (workflow/rules/03a_featureCounts.snake) for
# translational-efficiency analysis: only in-frame P-sites falling within
# the annotated CDS are counted (5'/3' UTR reads excluded), using riboWaltz's
# own cds_coverage() rather than whole-read exon overlap.
#
# Runs on the SAME bams/{sample}_Aligned.toTranscriptome.filtered.bam files
# produced by rules/02d_riboseq_transcriptome_filter.snake (multimappers2.py
# -k 1, grouped by read name -- a MAPQ-based filter does NOT work here, see
# that rule's own comment for why; -k 1 guarantees a surviving read has
# exactly one transcript-space alignment record, which is what makes summing
# per-transcript CDS counts up to gene level safe -- no read can appear
# under more than one transcript/gene after this filter).
#
# Reuses annotation_dt/psite_offset persisted by ribowaltz_psite.R instead
# of recomputing create_annotation()/psite() (the expensive, read-count-
# scaled steps there) a second time.

args <- commandArgs(trailingOnly = TRUE)

library("tidyverse")
library("riboWaltz")
library("rtracklayer")

bam_dir          <- args[1]                       # PROJ/bams (contains *_Aligned.toTranscriptome.filtered.bam)
samples          <- strsplit(args[2], " ")[[1]]   # sample names
gtf_file         <- args[3]                       # primary index's GTF -- for transcript->gene rollup
annotation_rds   <- args[4]                       # from ribowaltz_psite.R
psite_offset_rds <- args[5]                       # from ribowaltz_psite.R
index_name       <- args[6]                       # PRIMARY_INDEX, for output file naming (parity w/ featureCounts output)
outdir           <- args[7]                       # PROJ/counts

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

annotation_dt <- readRDS(annotation_rds)
psite_offset  <- readRDS(psite_offset_rds)

# riboWaltz's bamtolist() expects one bam per sample named <name>.bam in
# bamfolder -- same symlinking convention as ribowaltz_psite.R.
link_dir <- file.path(outdir, "psite_cds_bam_links")
dir.create(link_dir, showWarnings = FALSE)
for (s in samples) {
  src <- file.path(bam_dir, paste0(s, "_Aligned.toTranscriptome.filtered.bam"))
  dst <- file.path(link_dir, paste0(s, ".bam"))
  if (!file.exists(dst)) file.symlink(normalizePath(src), dst)
}

reads_list       <- bamtolist(bamfolder = link_dir, annotation = annotation_dt)
reads_psite_list <- psite_info(reads_list, psite_offset)

# One call across all samples at once -- cds_coverage() returns one row per
# transcript with an annotated CDS, one column per sample (transcript,
# length_cds, <sample1>, <sample2>, ...). in_frame = TRUE (default) restricts
# further to P-sites in the correct reading frame.
psite_cds <- cds_coverage(reads_psite_list, annotation_dt)

# ---- transcript -> gene rollup -----------------------------------------------
# annotation_dt (from create_annotation()) is transcript-only (no gene_id) --
# pull gene_id/gene_name/gene_biotype per transcript, and each gene's own
# Chr/Start/End/Strand, from the same GTF featureCounts already uses.
gtf <- rtracklayer::import(gtf_file)

tx_meta <- as.data.frame(gtf[gtf$type == "transcript"]) %>%
  dplyr::transmute(
    transcript_id,
    gene_id,
    gene_name = if ("gene_name" %in% names(.)) gene_name else gene_id,
    biotype   = if ("gene_biotype" %in% names(.)) gene_biotype else NA_character_
  )

gene_meta <- as.data.frame(gtf[gtf$type == "gene"]) %>%
  dplyr::transmute(
    gene_id,
    Chr    = as.character(seqnames),
    Start  = start,
    End    = end,
    Strand = as.character(strand)
  )

# riboWaltz's annotation_dt$transcript carries whatever identifier the
# transcriptome FASTA/BAM used (STAR's --quantMode TranscriptomeSAM
# reference names, derived from the same GTF) -- join on that directly
# against the GTF's own transcript_id attribute.
psite_cds_gene <- psite_cds %>%
  as.data.frame() %>%
  dplyr::rename(transcript_id = transcript) %>%
  dplyr::inner_join(tx_meta, by = "transcript_id") %>%
  dplyr::group_by(gene_id, gene_name, biotype) %>%
  dplyr::summarise(
    Length = max(length_cds, na.rm = TRUE),
    dplyr::across(all_of(samples), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::left_join(gene_meta, by = "gene_id")

# ---- per-sample output, matching featureCounts' column schema ---------------
# te_analysis.R's read_fc() only actually uses ID/gene_name/biotype/count --
# Chr/Start/End/Strand/Length are carried for schema parity (so this is a
# genuine drop-in alternative input for RiboSeq_TE.smk), not consumed.
for (s in samples) {
  out_tab <- psite_cds_gene %>%
    dplyr::transmute(
      ID         = gene_id,
      Chr, Start, End, Strand, Length,
      gene_name,
      biotype,
      count = .data[[s]]
    )
  write_tsv(out_tab, file.path(outdir, paste0(s, "_aligned_", index_name, "_psite_cds_featureCounts.tsv")))
}
