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
# res_dir <- "test/test_PRJNA753050"

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

# Gene regions
# named list with data.frames to intersect pauses with
regs <- plot_params$pause_regions %>%
  map(vroom, col_names = bed_cols)

# Subsample sample names
sub_regs <- set_names(plot_params$subsample_regions)

# Sample names
# the names in the samples list should correspond to the sampling group
sam_nms <- set_names(
  sample_df$sample,
  sample_df$file
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

# Treatments
n_sams <- length(sam_nms)

treat_reps <- sam_nms %>%
  str_split("_(?=[^_]*$)")

treats <- treat_reps %>%
  map_chr(pluck, 1) %>%
  unique()

n_treats <- length(treats)

# Replicates
reps <- treat_reps %>%
  map_chr(pluck, 2) %>%
  unique()

n_reps <- length(reps)

# Strings to match samples and create labels
clrs <- sample_df %>% distinct(clr) %>% unlist()
# Sample colors
clrs <- set_names(clrs, treats)


sub_pause_dir <- sub_sams %>%
  map(get_file_dir, "beds/pauses")

sub_meta_dir <- sub_sams %>%
  map(get_file_dir, "metaplot_beds")

p_dir <- sub_pause_dir$gene

list_dir <- file.path(plot_params$ref_dir, plot_params$list_dir)

# filtered lists
genes_pause <- load_genes(list_dir, plot_params$genes_pause, MIN_GENE_LEN)

genes <- gene_coords %>%
  mutate(pct = percent_rank(counts)) %>%
  dplyr::filter(counts > 0 & pct > plot_params$gene_min_pct) %>%
  dplyr::select(name, counts)

g_p <- genes_pause %>%
  semi_join(genes, by = "name")


# All pauses
pauses <- load_pause_beds(
  prfxs  = names(p_dir),
  sfxs   = str_c(pause_prfx, "pauses.bed.gz"),
  paths  = unname(p_dir),
  genes  = g_p,
  region = regs
) %>%
  format_sample_names(key_vec = sub_sam_lnms$gene)

# Create table of input parameters
nms <- bind_rows(
  expand_grid(treat = treats, region = names(regs)),
  expand_grid(rep   = reps,   region = names(regs))
)

pause_df <- pauses %>%
  nest(pauses = -c(sample, treat, rep, region))

# Create table with overlap for each comparison
pause_overlap <- nms %>%
  pmap_dfr(~ {
    filt_args <- list(...)
    filt_args <- filt_args[!is.na(filt_args)]
    
    # Filter data based on input arguments
    dat <- pause_df
    
    filt_args %>%
      iwalk(~ {
        dat <<- dat %>%
          dplyr::filter(!!sym(.y) == .x)
      })
    
    if (nrow(dat) < 2) {
      return(NULL)
    }
    
    # Calculate overlap
    dfs <- set_names(dat$pauses, dat$sample)
    res <- calc_pause_overlap(dfs)
    
    res <- dat %>%
      full_join(res, by = c(sample = "key")) %>%
      mutate(comp = str_c(names(dfs), collapse = "\n"))
    
    # Fill in NAs
    res <- res %>%
      mutate(across(everything(), ~ {
        if (n_distinct(.x, na.rm = TRUE) == 1) {
          replace_na(.x, unique(na.omit(.x)))
          
        } else {
          .x
        }
      }))
    
    res
  })

if (is_empty(pause_overlap)) {
  PLOT_OVERLAP <- FALSE
  n_comps <- 1
  
} else {
  PLOT_OVERLAP <- TRUE
  
  n_comps <- n_distinct(pause_overlap$comp)
  
  subtitle_text <- glue("<ul style='text-align: left;'>
                        <li>genes > {MIN_GENE_LEN} kb long and separated by > {extract_sep_info(plot_params$genes_pause)} kb
      are shown.</li>
                        <li>pauses were identified after downsampling entire gene region.</li>
                      </ul>")
  
  # Print pause stats
  my_table <- gt(pause_overlap %>%
    dplyr::filter(class == "unique") %>%
    distinct(sample, region, n_genes, n_pauses) %>%
    group_by(sample) %>%
    mutate(total_pauses = sum(n_pauses)) %>%
    ungroup()) %>% 
    tab_header(
      title = "The number of pause sites shared between samples and replicates is shown below for the TSS (TSS - +500 bp) and gene body (+500 bp - pAS) regions.",
      
      subtitle = html(subtitle_text)
    )
}

gtsave(my_table, snakemake@output[["overlap_table"]])

####
bars <- ggplot() +
  geom_blank()

if (PLOT_OVERLAP) {
  # Bar graph colors
  clr_df <- sample_df %>%
    separate(sample, sep = "_(?=[^_]*$)", into = c("treat", "rep"), remove = FALSE) %>%
    mutate(
      rep_i = as.numeric(str_extract(rep, "[0-9]+$")),
      clr   = clrs[treat],
      clr   = lighten(clr, 0.3 * (rep_i - 1))
    )
  
  bar_clrs <- set_names(clr_df$clr, clr_df$sample)
  bar_clrs <- c(shared = "grey85", bar_clrs)
  
  # Bar graph data
  plt_dat <- pause_overlap %>%
    mutate(
      sample = ifelse(class == "shared", class, sample),
      sample = fct_relevel(sample, names(bar_clrs)),
      comp   = fct_inorder(comp),
      region = fct_relevel(region, names(regs))
    ) %>%
    group_by(comp, region) %>%
    mutate(
      comp_reg = str_c(comp, "_", region),
      x_lab    = str_c(comma(sum(n)), "\npauses")
    ) %>%
    ungroup()
  
  # Fill labels for legend
  plt_dat <- plt_dat %>%
    group_by(comp) %>%
    mutate(
      fill_lab = case_when(
        n_distinct(treat, na.rm = TRUE) == 1 ~ rep,
        n_distinct(rep,   na.rm = TRUE) == 1 ~ treat
      ),
      fill_lab = if_else(
        is.na(fill_lab),
        as.character(sample),
        fill_lab
      )
    ) %>%
    ungroup()
  
  # X-axis labels
  x_labs <- distinct(plt_dat, comp_reg, x_lab)
  
  x_labs <- set_names(
    x_labs$x_lab,
    x_labs$comp_reg
  )
  
  # Split data.frame into list for plotting
  df_lst <- plt_dat %>%
    arrange(region) %>%
    mutate(comp_reg = fct_inorder(comp_reg)) %>%
    group_by(comp) %>%
    group_split()
  
  # Create bar graphs
  bars <- df_lst %>%
    map(~ {
      fill_labs <- set_names(
        unique(.x$fill_lab),
        unique(.x$sample)
      )
      
      plt_clrs <- bar_clrs[names(bar_clrs) %in% .x$sample]
      
      .x %>%
        mutate(
          sample = fct_drop(sample),
          sample = fct_relevel(sample, "shared", after = 1)
        ) %>%
        ggplot(aes(comp_reg, n, fill = sample)) +
        geom_col(position = "fill", width = 0.95) +
        scale_fill_manual(values = plt_clrs, labels = fill_labs) +
        scale_x_discrete(labels = x_labs) +
        facet_wrap(~ region, scales = "free") +
        
        ggtitle(unique(.x$comp)) +
        theme_info +
        theme(
          aspect.ratio     = 2.5,
          panel.border     = element_blank(),
          #panel.spacing    = unit(-10, "pt"),
          plot.title       = element_text(hjust = 0.5),
          legend.position  = "bottom",
          legend.title     = element_blank(),
          legend.direction = "vertical",
          strip.text       = element_text(vjust = 0),
          axis.title       = element_blank(),
          axis.line        = element_blank(),
          axis.ticks       = element_blank(),
          axis.text.x      = element_text(vjust = 4),
          axis.text.y      = element_blank()
        )
    })
  
  bars <- plot_grid(
    plotlist = bars,
    nrow     = 1,
    align    = "h",
    axis     = "tb"
  )+theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm")  # Top, right, bottom, left padding
  )
  
}

ggsave(snakemake@output[["overlap_bars"]], plot = bars, width = 12, height = 5)

####

# NET-seq signal is shown below for pause sites identified for TSS (TSS - +500 bp) and gene body (+500 bp - pAS) regions.
# The sequence preference for each group of sites is shown below each plot.
# The dotted line indicates the 3' end of the nascent RNA.
# 
# * genes >`r MIN_GENE_LEN` kb long and separated by >`r extract_sep_info(plot_params$genes_pause)` kb are shown (`r basename(plot_params$genes_pause)`)
# * pauses were identified after downsampling entire gene region
# * pauses separated by >30 bp are shown

# # Separated pauses
# # Do not format sample names until after merging with metaplot files
# sep_pauses <- load_pause_beds(
#   prfxs  = names(p_dir),
#   sfxs   = str_c(sep_prfx, "pauses.bed.gz"),
#   paths  = unname(p_dir),
#   genes  = g_p,
#   region = regs
# )
# 
# n_regs <- n_distinct(sep_pauses$region)
# Load and merge bed files
# sfxs <- c("sense" = str_c(pause_prfx, "pauses_meta_N.bed.gz"))
# grp <- "" #names(plot_grps$SAMPLES)
# win_cols <- c(bed_cols, "counts")
# win_cols[5] <- "win_id"
# 
# # out  <- here(plot_params$obj_dir, str_c(str_c(pause_prfx, "pauses_meta.tsv.gz")))
# out <- "/beevol/home/erickson/Mouse/NetSeq/230423_PRJNA753050/test/test_PRJNA753050/objects/test_PRJNA753050_200_strong_pauses_meta.tsv.gz"
# merge_p <- load_merge_wins(
#   prfxs         = names(p_dir),      # Sample names
#   sfxs          = sfxs,              # Suffix list with types as names
#   paths         = p_dir,             # Directories containing files
#   group         = grp,               # Group name
#   
#   file_out      = out,               # Path to output file
#   overwrite     = plot_params$overwrite,  # Overwrite output file if it exists
#   col_names     = win_cols,          # Column names for bed files
#   filter_unique = FALSE,             # Remove genes that are not shared between all samples
#   win_num       = 100,               # Total number of expected windows (including sense + anti)
#   win_min       = 41,
#   win_max       = 60,
#   ref_win       = NULL               # Reference window for adjusting window ids
# )
# 
# # Filter pauses for separated genes
# # add gene regions
# merge_p <- sep_pauses %>%
#   dplyr::select(name, gene_name, sample, region, n_lab) %>%
#   inner_join(merge_p, by = c("name", "sample"))
# 
# # Filter for overlapping genes
# # want all samples for each region to use the same genes
# # need to recalculate n values since filtering genes/pauses
# merge_p <- merge_p %>%
#   group_by(region) %>%
#   group_split() %>%
#   map_dfr(
#     merge_wins,
#     groups        = "sample",
#     ref_win       = NULL,
#     filter_unique = TRUE,
#     by            = "gene_name"
#   ) %>%
#   group_by(sample, region) %>%
#   mutate(
#     n_lab = comma(n_distinct(name)),
#     n_lab = str_c(n_lab, " pauses\n"),
#     n_lab = str_c(n_lab, comma(n_distinct(gene_name))),
#     n_lab = str_c(n_lab, " genes")
#   ) %>%
#   ungroup()
# 
# merge_clmns <- c("name", "gene_name", "sample", "region")
# 
# filt_pauses <- merge_p %>%
#   distinct(!!!syms(merge_clmns), n_lab)
# 
# # filter pauses so they match merge_p
# filt_pauses <- sep_pauses %>%
#   dplyr::select(-n_lab) %>%
#   inner_join(filt_pauses, by = merge_clmns) %>%
#   group_by(sample, region) %>%
#   mutate(
#     orig_lab = n_lab,
#     n_lab    = comma(n_distinct(name)),
#     n_lab    = str_c(n_lab, " pauses\n"),
#     n_lab    = str_c(n_lab, comma(n_distinct(gene_name))),
#     n_lab    = str_c(n_lab, " genes")
#   ) %>%
#   ungroup() %>%
#   format_sample_names(key_vec = sub_sam_lnms$gene)
# 
# stopifnot(all(filt_pauses$n_lab == filt_pauses$orig_lab))
# 
# filt_pauses <- filt_pauses %>%
#   dplyr::select(-orig_lab)
# 
# # Calculate mean signal
# mean_grps <- c("sample", "group", "type", "region")
# 
# mean_p <- merge_p %>%
#   calc_mean_pause_signal(
#     mean_grps,
#     key_vec = sub_sam_lnms$gene,
#     rel     = FALSE
#   )
# 
# rel_p <- merge_p %>%
#   calc_mean_pause_signal(
#     mean_grps,
#     key_vec = sub_sam_lnms$gene,
#     rel     = TRUE
#   )
# 
# # Create nested data.frame with data for metaplots and logos
# p_dat <- list(
#   "mean_counts" = mean_p,
#   "rel_counts"  = rel_p,
#   "pauses"      = filt_pauses
# )
# 
# sep_pause_df <- p_dat %>%
#   imap(~ {
#     .x %>%
#       nest(!!sym(.y) := -c(sample, region))
#   }) %>%
#   purrr::reduce(full_join, by = c("sample", "region")) %>%
#   separate(sample, sep = "_(?=[^_]*$)", into = c("treat", "rep"), remove = FALSE)
# 
# # Check for rows not shared between all data.frames
# if (!all(complete.cases(select_if(sep_pause_df, is.character)))) {
#   stop("Not all rows shared, malformed data.frame")
# }
# 
# # Check that n_lab matches for all data.frames
# n_chk <- sep_pause_df %>%
#   mutate(
#     across(
#       all_of(names(p_dat)),
#       ~ map_chr(.x, ~ unique(.x$n_lab))
#     )
#   ) %>%
#   pivot_longer(all_of(names(p_dat))) %>%
#   group_by(sample, region) %>%
#   dplyr::filter(length(unique(value)) > 1)
# 
# stopifnot(nrow(n_chk) == 0)
# 
# # Plot mean signal with logos
# y_lab    <- "RPKM"
# y_lim    <- range(mean_p$counts)
# y_lim[1] <- 0
# 
# # Genome files
# chrom_sizes <- read_genome(
#   file.path(plot_params$ref_dir, plot_params$chrom_sizes)
# )
# fa <- here(plot_params$ref_dir, plot_params$fasta)
# 
# logo_regs <- names(regs)
# logo_regs <- logo_regs[logo_regs %in% sep_pause_df$region]
# 
# mean_figs <- sep_pause_df %>%
#   mutate(
#     ttl = str_c(treat, rep, region, sep = " "),
#     fig = pmap(
#       list(mean_counts, pauses, ttl),
#       create_pause_logo_fig,
#       y_lab  = y_lab,
#       y_lim  = y_lim,
#       genome = chrom_sizes
#     ),
#     sample = fct_relevel(sample, sam_nms),
#     region = fct_relevel(region, logo_regs)
#   ) %>%
#   arrange(region, sample)
# 
# # Run this in separate chunk to exclude fasta messages
# plot_grid(
#   plotlist = mean_figs$fig,
#   ncol     = n_sams
# )
# 
# 
# 
# ### Relative signal
# 
# 
# y_lab    <- "relative signal"
# y_lim    <- range(rel_p$counts)
# y_lim[1] <- 0
# 
# rel_figs <- sep_pause_df %>%
#   mutate(
#     ttl = str_c(sample, " ", region),
#     fig = pmap(
#       list(rel_counts, pauses, ttl),
#       create_pause_logo_fig,
#       y_lab  = y_lab,
#       y_lim  = y_lim,
#       genome = chrom_sizes
#     ),
#     sample = fct_relevel(sample, sam_nms),
#     region = fct_relevel(region, logo_regs)
#   ) %>%
#   arrange(region, sample)
# 
# 
# 
# # Run this in separate chunk to exclude fasta messages
# plot_grid(
#   plotlist = rel_figs$fig,
#   ncol     = n_sams
# )

qpdf::pdf_combine(
  input = c(snakemake@output[["overlap_table"]], 
            snakemake@output[["overlap_bars"]]
            ),
  output = snakemake@output[["motifs"]]
)


