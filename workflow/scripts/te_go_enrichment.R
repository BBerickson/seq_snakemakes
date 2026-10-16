#!/usr/bin/env Rscript
# GO enrichment on TE-changed genes. Gated behind run_go_enrichment: false by
# default (see RiboSeq_TE_samples.yaml). clusterProfiler/AnnotationDbi/GO.db
# are genome-agnostic and live in the Ribo-seq container itself (confirmed
# present there) -- but the species-specific org.*.eg.db annotation package
# deliberately does NOT: it's kept on the server (bind-mounted, not baked
# into the image) in a shared R library directory, so adding support for a
# new genome is "install one more org.db package there", not a container
# rebuild. See ref/hg38_p14/go_annotation/README.txt for how this directory
# is set up and how to add another organism (e.g. mouse).

args <- commandArgs(trailingOnly = TRUE)

library("tidyverse")
library("clusterProfiler")

results_file  <- args[1]
orgdb_name    <- args[2]
outfile       <- args[3]
orgdb_lib_dir <- args[4]  # GO_ANNOTATION_LIB_DIR -- shared, non-container org.db location; "" to rely on .libPaths() as-is

if (nzchar(orgdb_lib_dir)) {
  .libPaths(c(orgdb_lib_dir, .libPaths()))
}

if (!requireNamespace(orgdb_name, quietly = TRUE)) {
  stop(sprintf(
    "GO enrichment requested but %s is not installed in GO_ANNOTATION_LIB_DIR (%s) or anywhere else on .libPaths() -- see ref/hg38_p14/go_annotation/README.txt to install it, or set run_go_enrichment: false.",
    orgdb_name, orgdb_lib_dir
  ))
}
library(orgdb_name, character.only = TRUE)

te <- read_tsv(results_file, show_col_types = FALSE)
sig_genes <- te %>%
  mutate(te_padj = coalesce(apvRvmPAdj, pvalue.adjust)) %>%
  dplyr::filter(!is.na(te_padj), te_padj < 0.05) %>%
  pull(ID)

universe <- te$ID

ego <- enrichGO(
  gene = sig_genes,
  universe = universe,
  OrgDb = get(orgdb_name),
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

write_tsv(as.data.frame(ego), outfile)
