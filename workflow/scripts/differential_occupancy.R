#!/usr/bin/env Rscript
# Differential ribosome occupancy: a plain two-group DESeq2 test on this
# Ribo-seq project's OWN featureCounts, no paired RNA-seq needed or used.
# This is NOT translational efficiency (see te_analysis.R for that, which
# needs a completed RNA-seq run) -- without RNA-seq there's no way to tell
# "more ribosomes on this transcript" apart from "more of this transcript
# exists", so a hit here could reflect transcriptional or translational
# change (or both). Still a real, standard thing to report when no paired
# RNA-seq run exists -- see RiboSeq_samples.yaml's CONTRASTS comment.
#
# Deliberately uses DESeq2's top-level DESeq() wrapper, not a lower-level
# call -- DESeq() has its own automatic parametric->local dispersion-fit
# fallback built in, which is exactly what xtail's direct low-level
# estimateDispersions() call lacks (confirmed the hard way: xtail crashed
# outright on a real 2-replicate-per-group contrast where gene-wise
# dispersions were too tightly clustered to curve-fit; see te_analysis.R's
# own comment). Still wrapped in tryCatch here for defense in depth, same
# pattern as te_analysis.R's xtail handling.

args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
})

ribo_files    <- strsplit(args[1], " ")[[1]]
group1        <- strsplit(args[2], " ")[[1]]
group2        <- strsplit(args[3], " ")[[1]]
contrast_name <- args[4]
filter_biotype<- args[5]
n_threads     <- as.integer(args[6])
results_tsv   <- args[7]
ma_png        <- args[8]
volcano_png   <- args[9]

samples <- c(group1, group2)
condition <- factor(c(rep("group1", length(group1)), rep("group2", length(group2))),
                     levels = c("group1", "group2"))

# featureCounts column convention shared with te_analysis.R/results.R/
# results_riboseq.R: Geneid, Chr, Start, End, Strand, Length, gene_name,
# biotype, count
read_fc <- function(file, sample) {
  col_names <- c("ID", "Chr", "Start", "End", "Strand", "Length", "gene_name", "biotype", "count")
  dat <- read_tsv(file, comment = "#", show_col_types = FALSE)
  names(dat) <- col_names
  dat %>% dplyr::select(ID, gene_name, biotype, !!sample := count)
}

# Snakemake's expand() already produced ribo_files in exactly samples'
# order (group1 then group2) -- no re-matching by filename needed.
tabs <- map2(ribo_files, samples, ~read_fc(.x, .y))
ribo_tab <- purrr::reduce(tabs, full_join, by = c("ID", "gene_name", "biotype"))
if (nzchar(filter_biotype)) {
  ribo_tab <- ribo_tab %>% dplyr::filter(biotype == filter_biotype)
}

gene_meta <- ribo_tab %>% dplyr::select(ID, gene_name, biotype)
ribo_mat <- ribo_tab %>% dplyr::select(all_of(samples)) %>% as.matrix()
rownames(ribo_mat) <- gene_meta$ID

# drop genes with zero counts everywhere -- DESeq2 chokes on all-zero rows
keep <- rowSums(ribo_mat) > 0
ribo_mat <- ribo_mat[keep, , drop = FALSE]
gene_meta <- gene_meta[keep, ]

res_table <- tryCatch({
  dds <- DESeqDataSetFromMatrix(
    countData = ribo_mat,
    colData = data.frame(condition = condition, row.names = samples),
    design = ~condition
  )
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "group2", "group1"), parallel = n_threads > 1)
  as.data.frame(res) %>% rownames_to_column("ID")
}, error = function(e) {
  warning(sprintf(
    "Differential occupancy test failed for %s: %s. Writing gene list with no effect/significance columns.",
    contrast_name, conditionMessage(e)
  ))
  gene_meta %>%
    dplyr::select(ID) %>%
    dplyr::mutate(baseMean = NA_real_, log2FoldChange = NA_real_, pvalue = NA_real_, padj = NA_real_)
})

combined <- gene_meta %>%
  left_join(res_table, by = "ID") %>%
  mutate(sig = !is.na(padj) & padj < 0.05)

write_tsv(combined, results_tsv)

# ---- plots -------------------------------------------------------------------
ma_plot <- ggplot(combined, aes(x = log2(baseMean + 1), y = log2FoldChange, color = sig)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey60", "firebrick")) +
  labs(x = "log2(mean footprint count + 1)", y = "Ribosome occupancy log2 fold change",
       title = paste0(contrast_name, ": differential occupancy MA plot"), color = "padj < 0.05") +
  theme_minimal()
ggsave(ma_png, ma_plot, width = 7, height = 5)

volcano_plot <- ggplot(combined, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("grey60", "firebrick")) +
  labs(x = "Ribosome occupancy log2 fold change", y = "-log10(adjusted p-value)",
       title = paste0(contrast_name, ": differential occupancy volcano plot"), color = "padj < 0.05") +
  theme_minimal()
ggsave(volcano_png, volcano_plot, width = 7, height = 5)
