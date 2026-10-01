#!/usr/bin/env Rscript
# GO enrichment on TE-changed genes. Gated behind run_go_enrichment: false by
# default (see RiboSeq_TE_samples.yaml) since clusterProfiler + an org.*.eg.db
# annotation package are NOT in the Ribo-seq container's
# install_packages_riboseq.R as pasted this session -- add them there first.

args <- commandArgs(trailingOnly = TRUE)

library("tidyverse")
library("clusterProfiler")

results_file <- args[1]
orgdb_name   <- args[2]
outfile      <- args[3]

if (!requireNamespace(orgdb_name, quietly = TRUE)) {
  stop(sprintf(
    "GO enrichment requested but %s is not installed -- add it (and clusterProfiler) to install_packages_riboseq.R, or set run_go_enrichment: false.",
    orgdb_name
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
