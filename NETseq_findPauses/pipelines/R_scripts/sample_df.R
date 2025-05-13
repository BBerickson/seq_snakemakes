#!/usr/bin/env Rscript

library(tidyverse)
library(yaml)

# setwd("~/Mouse/NetSeq/230423_PRJNA753050/test")
sam_yam <- snakemake@input[["sam_yam"]]
# sam_yam <- "samples.yaml"
myfuns <- snakemake@params[["myfuns"]]
# myfuns <- "pipelines/R_scripts/funs.R"
setup <- snakemake@params[["setup"]]
# setup <- "pipelines/R_scripts/setup.R"

source(setup)
source(myfuns)

plot_grps <- yaml::read_yaml(sam_yam)
# plot_grps <- yaml::read_yaml("samples.yaml")

# Set missing sample names
# set default names if not provided
plot_grps$SAMPLES <- plot_grps$SAMPLES %>%
  map(~ {
    .x %>%
      map(~ {
        if (is.null(names(.x))) names(.x) <- str_c("sample-", seq_along(.x))
        
        as.list(.x)
      })
  })

# color setup
sections   <- names(plot_grps$SAMPLES)
plot_grps$COLORS <- set_names(list(unique(plot_grps$COLORS)),sections)
sec_colors <- names(plot_grps$COLORS)

missing_colors <- sections[!sections %in% sec_colors]

plot_grps$COLORS[missing_colors] <- plot_grps$SAMPLES[missing_colors] %>%
  map(~ {
    n_clrs <- .x %>%
      map(names) %>%
      unlist(use.names = FALSE) %>%
      n_distinct()
    
    c("black", theme_colors[seq_len(n_clrs - 1)])
  })

# convert RGB to hex and makes sure only 2 colors
for(grp in names(plot_grps$SAMPLES)){
  plot_grps$COLORS[[grp]] <- plot_grps$COLORS[[grp]][1:2] %>%
    map(~ {
      RgbToHex(.x)
    }) %>% unlist()
}

# Set replicate names for plots
# set default replicate IDs if not provided _r#
plot_grps$SAMPLES <- plot_grps$SAMPLES %>%
  map(~ {
    ids <- str_c("_r", seq_along(.x))
    
    map2(.x, ids, ~ {
      names(.x) <- str_c(names(.x), .y)
      
      as.list(.x)
    })
  })

# Create data.frame with sample info
sample_df <- plot_grps$SAMPLES %>%
  imap_dfr(~ {
    plt_grp <- .y
    clrs    <- plot_grps$COLORS[[plt_grp]]
    
    imap_dfr(.x, ~ {
      grp <- .y
      sam_clrs <- set_names(clrs, names(.x))
      
      imap_dfr(.x, ~ {
        sam <- .y
        
        map_dfr(.x, ~ {
          file <- .x
          
          tibble(
            sampling_grp = grp,
            plot_grp     = plt_grp,
            sample       = sam,
            file         = file,
            clr          = sam_clrs[[sam]]
          )
        })
      })
    })
  })


sample_df <- sample_df %>%
  mutate(sample_parts = strsplit(sample, "_|-")) %>%
  unnest_wider(sample_parts, names_sep = "_") %>%
  rename_with(~ "rep", last(starts_with("sample_parts_"))) %>%
  group_by(clr) %>%
  mutate(across(starts_with("sample_parts_"), ~ ifelse(all(. == .[1]), ., ""))) %>%
  mutate(across(starts_with("sample_parts_"), ~ na_if(., ""))) %>%
  unite("sample", starts_with("sample_parts_"), sep = "_", na.rm = TRUE) %>%
  mutate(sample = ifelse(is.na(sample) | sample == "", "sample", sample)) %>%
  ungroup() %>% 
  mutate(sam=sample,sample=paste0(sample,"_",rep))

#### save file
write_tsv(sample_df, snakemake@output[[1]])
