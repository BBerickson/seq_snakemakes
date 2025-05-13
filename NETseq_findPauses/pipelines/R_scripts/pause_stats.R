#!/usr/bin/env Rscript

library(tidyverse)
library(yaml)

# setwd("~/Mouse/NetSeq/230423_PRJNA753050/test")

plot_yam <- snakemake@input[["plot_yam"]]
# plot_yam <- "pipelines/plots.yaml"
sample_df <- snakemake@input[["sample_df"]]
# sample_df <- "test_PRJNA753050/objects/sample_df.tsv"
gene_coords <- snakemake@input[["gene_coords"]]
# gene_coords <- "test_PRJNA753050/objects/metaplot_genes.tsv.gz"
genome_params <- snakemake@params[["genome_params"]]
# genome_params <- "pipelines/ref/mm39_FP.yaml"
myfuns <- snakemake@params[["myfuns"]]
# myfuns <- "pipelines/R_scripts/funs.R"
setup <- snakemake@params[["setup"]]
# setup <- "pipelines/R_scripts/setup.R"
res_dir <- snakemake@params[["res_dir"]]
# res_dir <- "test_PRJNA753050"

# load functions and files
source(setup)
source(myfuns)


plot_params <- yaml::read_yaml(plot_yam)
genome_params <- yaml::read_yaml(genome_params)
plot_params   <- append(plot_params, genome_params)
sample_df <- read_tsv(sample_df)
gene_coords <- read_tsv(gene_coords)

# Pause parameters
win   <- plot_params$pause_win
stren <- plot_params$pause_strength
pause_prfx <- str_c("_", win, stren)
sep_prfx   <- str_c("_", win, "_sep", stren)

# filtered lists
list_dir <- file.path(plot_params$ref_dir, plot_params$list_dir)
genes_pause <- load_genes(list_dir, plot_params$genes_pause, MIN_GENE_LEN)
genes <- gene_coords %>%
  mutate(pct = percent_rank(counts)) %>%
  dplyr::filter(counts > 0 & pct > plot_params$gene_min_pct) %>%
  dplyr::select(name, counts)
g_p <- genes_pause %>%
  semi_join(genes, by = "name")

# Subsample sample names
sub_regs <- set_names(plot_params$subsample_regions)

# Load gene region pause bed files
regs <- c(
  "5-TSS-100bp"   = "TSS-100bp",
  "5-100bp-300bp" = "100bp-300bp",
  "5-300bp-500bp" = "300bp-500bp",
  "5-500bp-1kb"   = "500bp-1kb",
  "tss"           = "tss",
  "body"          = "body",
  "gene"          = "gene"
)

sam_lnms <- set_names(
  sample_df$sample,
  str_c(sample_df$file, "-", sample_df$sampling_grp)
)

sams <- names(sam_lnms)

sub_sams <- sub_regs %>%
  map(~ str_c(sams, "_", .x))

sub_sam_lnms <- sub_sams %>%
  map(~ set_names(sample_df$sample, .x))

gene_beds <- sub_regs %>%
  map(~ {
    crossing(
      sam = sub_sams[[.x]],
      reg = names(regs)
    )
  })

gene_beds <- gene_beds %>%
  imap(~ {
    .x %>%
      pmap_dfr(load_region_beds, genes = g_p, prfx = pause_prfx) %>%
      mutate(len = end - start) %>%
      format_sample_names(key_vec = sub_sam_lnms[[.y]]) %>%
      dplyr::select(-c(chrom, start, end, score, strand)) %>%
      mutate(sub_reg = .y) %>%
      relocate(type, counts, .after = last_col())
  })

# Combine data.frames for TSS and body subsampling regions
# * TSS windows use gene subsampling data
# * TSS region uses TSS subsampling data
# * body region uses body subsampling data
gene_df <- gene_beds$gene %>%
  dplyr::filter(!region %in% sub_regs)

tss_df <- gene_beds$tss %>%
  dplyr::filter(region == "tss")

body_df <- gene_beds$body %>%
  dplyr::filter(region == "body")

gene_beds <- bind_rows(gene_df, tss_df, body_df)

# Check data.frame for obvious issues
pause_stats <- gene_beds %>%
  pivot_wider(names_from = "type", values_from = "counts")

chk1 <- pause_stats %>%
  dplyr::filter(pauses > 0 & (`p-reads` == 0 | NET == 0))  # regions with pauses but no pause reads

chk2 <- pause_stats %>%
  dplyr::filter(`p-reads` > NET)                           # region with more pause reads than total reads

chk3 <- pause_stats %>%                                    # bad subsample region
  dplyr::filter(
    (sub_reg == "tss" & sub_reg != region) |
      (sub_reg == "body" & sub_reg != region)
  )

stopifnot(
  nrow(chk1) == 0,
  nrow(chk2) == 0,
  nrow(chk3) == 0
)

# Calculate pausing stats
# remove regions with 0 mNET-seq reads
pause_stats <- pause_stats %>%
  dplyr::filter(NET > 0) %>%
  calc_pause_stats() %>%
  dplyr::select(-len, -sub_reg)

pause_stats %>%
  write_tsv(snakemake@output[[1]])

