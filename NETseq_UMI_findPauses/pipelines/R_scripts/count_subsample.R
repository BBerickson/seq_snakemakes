#!/usr/bin/env Rscript

library("tidyverse")

samp <- list()
for(i in seq_along(snakemake@input)) {
  samp[[i]] <- read_delim(snakemake@input[[i]],col_names = c("name","Filtered","value"),delim = " ") 
  }

out <- bind_rows(samp) %>% dplyr::mutate(value=min(value)/value, group=snakemake@params[["group"]]) %>% dplyr::select(name, group, value)

write_tsv(out, snakemake@output[[1]],col_names = F)
