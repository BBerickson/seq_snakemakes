#!/usr/bin/env Rscript

library(tidyverse)
library(yaml)

# setwd("~/Mouse/NetSeq/230423_PRJNA753050/test")

pause_stats <- snakemake@input[["pause_stats"]]
# pause_stats <- "test_PRJNA753050/objects/pause_stats.tsv.gz"
plot_yam <- snakemake@input[["plot_yam"]]
# plot_yam <- "pipelines/plots.yaml"
sample_df <- snakemake@input[["sample_df"]]
# sample_df <- "test_PRJNA753050/objects/sample_df.tsv"
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
pause_stats <- read_tsv(pause_stats)

# Sample names
# the names in the samples list should correspond to the sampling group
sam_nms <- set_names(
  sample_df$sample,
  sample_df$file
)

# Treatments
treat_reps <- sam_nms %>%
  str_split("_(?=[^_]*$)")

treats <- treat_reps %>%
  map_chr(pluck, 1) %>%
  unique()
n_treats <- length(treats)

common_reps <- sample_df %>%
  group_by(rep) %>%
  dplyr::filter(all(treats %in% sam)) %>%
  pull(rep) %>%
  unique()

n_common_reps <- length(common_reps)

# Treatment groups to use when calculating fold changes
fc_treats <- c(treats[n_treats - 1], last(treats))

if (identical(n_treats, 1L)) {
  fc_treats <- rep(fc_treats, 2)
}

# Strings to match samples and create labels
clrs <- sample_df %>% distinct(clr) %>% unlist()
clrs <- set_names(clrs, treats)

# Minimum number of pauses to use for filtering TSS boxplots and scatter plots
# First number is for boxplots, second number is for scatter plots
PAUSE_LIM <- plot_params$pause_lims

# Cutoffs for filtering TSS regions and body region for boxplots
TSS_LIM     <- plot_params$tss_lim      # raw number of reads
BODY_LIM    <- plot_params$body_lim     # reads / kb
ONLY_SHARED <- plot_params$only_shared  # include only shared genes

# gene and pause regions
regs <- c(
  "5-TSS-100bp"   = "TSS-100bp",
  "5-100bp-300bp" = "100bp-300bp",
  "5-300bp-500bp" = "300bp-500bp",
  "5-500bp-1kb"   = "500bp-1kb",
  "tss"           = "tss",
  "body"          = "body",
  "gene"          = "gene"
)

metrics <- unique(pause_stats$type)

# Filter regions based on signal cutoffs
box_dat_0 <- pause_stats %>%
  get_box_data(
    p_lim       = PAUSE_LIM[1],
    tss_lim     = TSS_LIM,
    body_lim    = BODY_LIM,
    only_shared = ONLY_SHARED
  )

# If ONLY_SHARED, do not require all regions to meet pause cutoff for box_dat_1
# * ONLY_SHARED should never be applied to pause cutoffs
# * ONLY_SHARED should be applied to other cutoffs
# * need to check for overlapping genes for samples since get_box_data does not
#   filter based on NAs
if (ONLY_SHARED) {
  box_dat_1 <- pause_stats %>%
    get_box_data(
      p_lim       = 0,        # set to 0 since we do not want to apply
      tss_lim     = TSS_LIM,  # ONLY_SHARED to pause cutoffs
      body_lim    = BODY_LIM,
      only_shared = ONLY_SHARED
    )
  
  # Still need to ensure that plotted genes are shared between samples
  box_dat_1 <- box_dat_1 %>%
    pivot_wider(names_from  = "type", values_from = "counts") %>%
    dplyr::filter(pauses >= PAUSE_LIM[2]) %>%
    pivot_longer(all_of(metrics), names_to = "type", values_to = "counts") %>%
    
    pivot_wider(names_from = "treat", values_from = "counts") %>%
    dplyr::filter(if_all(all_of(treats), ~ !is.na(.x))) %>%
    pivot_longer(all_of(treats), names_to = "treat", values_to = "counts")
  
} else {
  box_dat_1 <- pause_stats %>%
    get_box_data(
      p_lim       = PAUSE_LIM[2],
      tss_lim     = TSS_LIM,
      body_lim    = BODY_LIM,
      only_shared = ONLY_SHARED
    )
}

# For strength vs density scatter plots do not require shared regions
box_dat_sctr <- pause_stats %>%
  get_box_data(
    p_lim       = PAUSE_LIM[2],
    tss_lim     = TSS_LIM,
    body_lim    = BODY_LIM,
    only_shared = FALSE
  )

# Regions to exclude from boxplots
boxes_exclude_regs <- "tss"

metric <- c("fraction pause reads" = "p-reads_NET")

# Only use reps that are present for all samples
p <- common_reps %>%
  map(~ {
    box_dat_0 %>%
      dplyr::filter(
        !region %in% boxes_exclude_regs,
        rep == .x
      ) %>%
      create_stat_boxes(metric, reg_labels = regs) +
      scale_y_continuous(
        breaks = seq(0, 1, 0.5),
        expand = expansion(mult = c(0.05, 0.2))
      ) +
      labs(title = .x)
  }) %>%
  plot_grid(plotlist = ., nrow = 1)

ggsave(snakemake@output[["fpr"]],p, width = 12, height = 5)

# The fraction of reads aligning to pauses is shown for TSS and gene body regions. p-values were calculated using the Wilcoxon rank sum test.
# 
# * genes >`r MIN_GENE_LEN` kb long and separated by >`r extract_sep_info(plot_params$genes_pause)` kb
# are shown (`r basename(plot_params$genes_pause)`)
# * TSS region downsampled for TSS boxplot
# * body region downsampled for gene body boxplot
# * TSS regions >=`r TSS_LIM` reads
# * body region >=`r BODY_LIM` read/kb
# * TSS and body regions >= `r PAUSE_LIM[2]` pause

# Data for boxplots
metric <- c("fraction pause reads" = "p-reads_NET")

box_dat <- pause_stats %>%
  filter(region %in% c("tss", "body")) %>%
  get_box_data(
    p_lim       = PAUSE_LIM[2],
    tss_lim     = TSS_LIM,
    body_lim    = BODY_LIM,
    only_shared = ONLY_SHARED
  )

# Only use reps that are present for all samples
p <- common_reps %>%
  map(~ {
    box_dat %>%
      dplyr::filter(rep == .x) %>%
      create_stat_boxes(metric, reg_labels = regs) +
      scale_y_continuous(
        breaks = seq(0, 1, 0.5),
        expand = expansion(mult = c(0.05, 0.2))
      ) +
      labs(title = .x) +
      theme(aspect.ratio = 1.5)
  }) %>%
  plot_grid(plotlist = ., nrow = 1)

ggsave(snakemake@output[["fprtb"]],p, width = 12, height = 5)

# The average number of reads aligning to each pause site is shown for regions described above that have at least one pause detected in the region.
# 
# * genes >`r MIN_GENE_LEN` kb long and separated by >`r extract_sep_info(plot_params$genes_pause)` kb
# are shown (`r basename(plot_params$genes_pause)`)
# * entire gene region downsampled for TSS regions
# * body region downsampled for gene body
# * TSS regions >=`r TSS_LIM` reads
# * body region >=`r BODY_LIM` read/kb
# * all regions >= `r PAUSE_LIM[2]` pause (there must be at least 1 pause to calculate reads / pause)
# 

metric <- c("reads / pause" = "p-reads_pauses")

# By default for p-reads_pauses regions will have >0 pauses due to NAs during
# calculation
p <- common_reps %>%
  map(~ {
    box_dat_1 %>%
      dplyr::filter(
        !region %in% boxes_exclude_regs,
        rep == .x
      ) %>%
      create_stat_boxes(metric, y_trans = "log10", reg_labels = regs) +
      labs(title = .x)
  }) %>%
  plot_grid(plotlist = ., nrow = 1)

ggsave(snakemake@output[["frp"]],p, width = 12, height = 5)
