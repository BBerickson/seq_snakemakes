#!/usr/bin/env Rscript
# Translational-efficiency (TE) analysis joining a completed Ribo-seq run's
# and a completed RNA-seq run's featureCounts outputs. Primary: anota2seq
# (separates translational/transcriptional/buffering regulation). Secondary:
# xtail (independent statistical approach, cross-check). Follows this repo's
# commandArgs + tidyverse convention (see results.R/fragment_results.R), with
# base R/S4 object handling where the Bioconductor packages require it.
#
# NOTE: anota2seq's/xtail's exact accessor function names below were verified
# against their public documentation (not executed against a live R session
# in this environment, since no R/Bioconductor environment is available
# here) -- treat this as a solid first draft to debug against your real
# container, flagged clearly at the one spot (anota2seqGetOutput) most likely
# to need a small adjustment.

args <- commandArgs(trailingOnly = TRUE)

library("tidyverse")
library("anota2seq")
library("xtail")

ribo_files    <- strsplit(args[1], " ")[[1]]
rna_files     <- strsplit(args[2], " ")[[1]]
group1        <- strsplit(args[3], " ")[[1]]
group2        <- strsplit(args[4], " ")[[1]]
contrast_name <- args[5]
filter_biotype<- args[6]
n_threads     <- as.integer(args[7])
outdir        <- args[8]

samples <- c(group1, group2)
condition <- factor(c(rep("group1", length(group1)), rep("group2", length(group2))),
                     levels = c("group1", "group2"))

# featureCounts column convention shared with results.R/results_riboseq.R:
# Geneid, Chr, Start, End, Strand, Length, gene_name, biotype, count
read_fc <- function(file, sample) {
  # same read-then-rename pattern as results.R/results_riboseq.R: comment="#"
  # skips featureCounts' leading "# Program:..." line, the real header row
  # (Geneid/Chr/.../gene_biotype/<bam path>) becomes the column names, which
  # are then overridden positionally.
  col_names <- c("ID", "Chr", "Start", "End", "Strand", "Length", "gene_name", "biotype", "count")
  dat <- read_tsv(file, comment = "#", show_col_types = FALSE)
  names(dat) <- col_names
  dat %>% dplyr::select(ID, gene_name, biotype, !!sample := count)
}

build_matrix <- function(files, samples) {
  tabs <- map2(files, samples, ~read_fc(.x, .y))
  merged <- reduce(tabs, full_join, by = c("ID", "gene_name", "biotype"))
  if (nzchar(filter_biotype)) {
    merged <- merged %>% dplyr::filter(biotype == filter_biotype)
  }
  merged
}

# Snakemake's expand() already produced ribo_files/rna_files in exactly
# samples' order (group1 then group2) -- no re-matching by filename needed.
ribo_tab <- build_matrix(ribo_files, samples)
rna_tab  <- build_matrix(rna_files, samples)

gene_meta <- ribo_tab %>% dplyr::select(ID, gene_name, biotype)
ribo_mat <- ribo_tab %>% dplyr::select(all_of(samples)) %>% as.matrix()
rna_mat  <- rna_tab  %>% dplyr::select(all_of(samples)) %>% as.matrix()
rownames(ribo_mat) <- gene_meta$ID
rownames(rna_mat)  <- gene_meta$ID

# drop genes with zero counts everywhere in either matrix -- both tools
# choke on all-zero rows
keep <- rowSums(ribo_mat) > 0 & rowSums(rna_mat) > 0
ribo_mat <- ribo_mat[keep, , drop = FALSE]
rna_mat  <- rna_mat[keep, , drop = FALSE]
gene_meta <- gene_meta[keep, ]

# ---- anota2seq (primary) ----------------------------------------------------
a2s <- anota2seqDataSetFromMatrix(
  dataP = ribo_mat,
  dataT = rna_mat,
  phenoVec = as.character(condition),
  dataType = "RNAseq",
  normalize = TRUE
)
a2s <- anota2seqRun(a2s, performQC = FALSE, performROT = FALSE)

# anota2seqGetOutput is the documented generic accessor for pulling a full
# results table out of the run object -- if your anota2seq version exposes
# this slightly differently, this is the one line to adjust.
a2s_translation <- anota2seqGetOutput(a2s, output = "full", selContrast = 1, analysis = "translation")
a2s_translation <- as.data.frame(a2s_translation) %>% rownames_to_column("ID")

# ---- xtail (secondary / complementary) --------------------------------------
xt <- xtail(rna_mat, ribo_mat, as.character(condition), baseLevel = "group1", threads = n_threads)
xt_table <- resultsTable(xt, log2FCs = TRUE) %>% rownames_to_column("ID")
# resultsTable(log2FCs=TRUE)'s exact column names (expected: mRNA/RPF log2FC
# columns alongside the TE log2FC/pvalue.adjust columns used below) weren't
# verified against a live run -- print them once so a mismatch is obvious
# immediately rather than silently producing blank plot axes.
message("xtail resultsTable columns: ", paste(names(xt_table), collapse = ", "))

# ---- combine -----------------------------------------------------------------
combined <- gene_meta %>%
  left_join(a2s_translation, by = "ID") %>%
  left_join(xt_table, by = "ID", suffix = c("_anota2seq", "_xtail"))

write_tsv(combined, file.path(outdir, paste0(contrast_name, "_TE_results.tsv")))

# ---- plots -------------------------------------------------------------------
# column names are best-effort based on anota2seq/xtail's documented output
# shape (apvEff = translation effect size, apvRvmPAdj = adjusted p-value for
# anota2seq; log2FC_TE_v2/pvalue.adjust for xtail) -- adjust if your installed
# versions name these differently.
plot_dat <- combined %>%
  mutate(
    te_effect = coalesce(apvEff, log2FC_TE_v2),
    te_padj   = coalesce(apvRvmPAdj, pvalue.adjust),
    sig       = !is.na(te_padj) & te_padj < 0.05
  )

ma_plot <- ggplot(plot_dat, aes(x = log2(rowMeans(rna_mat) + 1), y = te_effect, color = sig)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey60", "firebrick")) +
  labs(x = "log2(mean RNA-seq count + 1)", y = "TE log2 fold change",
       title = paste0(contrast_name, ": TE MA plot"), color = "padj < 0.05") +
  theme_minimal()
ggsave(file.path(outdir, paste0(contrast_name, "_MA_plot.png")), ma_plot, width = 7, height = 5)

volcano_plot <- ggplot(plot_dat, aes(x = te_effect, y = -log10(te_padj), color = sig)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey60", "firebrick")) +
  labs(x = "TE log2 fold change", y = "-log10(adjusted p-value)",
       title = paste0(contrast_name, ": TE volcano plot"), color = "padj < 0.05") +
  theme_minimal()
ggsave(file.path(outdir, paste0(contrast_name, "_volcano_plot.png")), volcano_plot, width = 7, height = 5)

quad_plot <- ggplot(plot_dat, aes(x = log2FC_mRNA, y = log2FC_RPF, color = sig)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("grey60", "firebrick")) +
  labs(x = "RNA-seq log2 fold change", y = "Ribo-seq log2 fold change",
       title = paste0(contrast_name, ": 4-quadrant plot"), color = "padj < 0.05") +
  theme_minimal()
ggsave(file.path(outdir, paste0(contrast_name, "_4quadrant_plot.png")), quad_plot, width = 6, height = 6)
