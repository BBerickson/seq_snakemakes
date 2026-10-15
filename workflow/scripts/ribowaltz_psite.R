#!/usr/bin/env Rscript
# P-site offset calibration, periodicity, metagene, and footprint length-
# distribution plots via riboWaltz. Follows this repo's established R-script
# convention (commandArgs + tidyverse), not the orphaned snakemake@input-
# object style used only by the unused featurecount_PCA.R.
#
# One invocation per project (not per-sample) -- riboWaltz's own API loads
# every sample's transcriptome BAM together as one list, so P-site offsets
# and periodicity are computed jointly across all samples in this run.

args <- commandArgs(trailingOnly = TRUE)

library("tidyverse")
library("riboWaltz")

bam_dir  <- args[1]                       # PROJ/bams (contains *_Aligned.toTranscriptome.out.bam)
samples  <- strsplit(args[2], " ")[[1]]   # sample names (bam basenames, minus the Aligned.toTranscriptome.out.bam suffix)
gtf_file <- args[3]                       # primary index's GTF
fa_file  <- args[4]                       # primary index's FA_SAMPLE (unused by riboWaltz itself here, kept for parity with other scripts' arg lists)
outdir   <- args[5]                       # PROJ/report/riboWaltz -- per-sample PNGs go here
outfile  <- args[6]                       # PROJ/stats/PROJ_ribowaltz_psite.tsv -- one summary row per sample

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# riboWaltz's bamtolist() expects one bam per sample named <name>.bam in
# bamfolder -- symlink the fixed STAR transcriptome-BAM names into riboWaltz's
# expected <sample>.bam naming.
link_dir <- file.path(outdir, "bam_links")
dir.create(link_dir, showWarnings = FALSE)
for (s in samples) {
  src <- file.path(bam_dir, paste0(s, "_Aligned.toTranscriptome.out.bam"))
  dst <- file.path(link_dir, paste0(s, ".bam"))
  if (!file.exists(dst)) file.symlink(normalizePath(src), dst)
}

annotation_dt <- create_annotation(gtfpath = gtf_file)
reads_list    <- bamtolist(bamfolder = link_dir, annotation = annotation_dt)

psite_offset     <- psite(reads_list, flanking = 6, extremity = "auto")
reads_psite_list <- psite_info(reads_list, psite_offset)

# one summary row per sample: P-site offset, in-frame %, dominant footprint length
summary_rows <- NULL
for (s in samples) {
  # frame_psite()'s default plot_style ("split") returns "count_dt" (raw
  # dt key) and one "plot_<sample>" entry per sample (not a flat "dt"/
  # "plot" key) -- same for rlength_distr()/metaprofile_psite() below.
  frames <- frame_psite(reads_psite_list, annotation_dt, sample = s, region = "cds")
  frame_tab <- frames[["count_dt"]]
  in_frame_pct <- frame_tab %>%
    dplyr::filter(frame == 0) %>%
    dplyr::summarise(pct = sum(scaled_count, na.rm = TRUE)) %>%
    dplyr::pull(pct)

  length_dist <- rlength_distr(reads_list, sample = s)
  length_tab  <- length_dist[["count_dt"]]
  dominant_len <- length_tab %>%
    dplyr::slice_max(count, n = 1, with_ties = FALSE) %>%
    dplyr::pull(length)

  # psite_offset has one row per (sample, read length); offset_from_5 is
  # frequently NA for rare lengths (not enough reads to estimate it), so
  # grabbing row 1 (sorted by length ascending) was returning NA for every
  # sample. corrected_offset_from_5 is riboWaltz's own backfilled/adjusted
  # value (never NA, and what psite_info() actually assigns to reads), and
  # tying it to the dominant (most common) footprint length -- rather than
  # whatever length happens to sort first -- gives a meaningful, stable
  # representative offset per sample.
  offset_row <- psite_offset %>%
    dplyr::filter(sample == !!s, length == dominant_len[1])

  summary_rows <- dplyr::bind_rows(
    summary_rows,
    tibble::tibble(
      sample = s,
      psite_offset = if (nrow(offset_row) > 0) offset_row$corrected_offset_from_5[1] else NA_real_,
      in_frame_pct = if (length(in_frame_pct) > 0) in_frame_pct else NA_real_,
      dominant_footprint_length = if (length(dominant_len) > 0) dominant_len else NA_real_
    )
  )

  # per-sample plots
  periodicity <- frames[[paste0("plot_", s)]]
  ggsave(file.path(outdir, paste0(s, "_periodicity.png")), periodicity, width = 7, height = 5)

  metaprofile <- metaprofile_psite(reads_psite_list, annotation_dt, sample = s)
  ggsave(file.path(outdir, paste0(s, "_metaprofile.png")), metaprofile[[paste0("plot_", s)]], width = 7, height = 5)

  ggsave(file.path(outdir, paste0(s, "_length_dist.png")), length_dist[[paste0("plot_", s)]], width = 7, height = 5)
}

write_tsv(summary_rows, outfile)
