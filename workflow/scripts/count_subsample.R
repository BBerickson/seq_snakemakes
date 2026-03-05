#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library("tidyverse")

samp_file <- strsplit(args[1], " ")[[1]]
mygroup <- args[2]
outfile <- args[3]

samp <- list()
for(i in seq_along(samp_file)) {
  samp[[i]] <- read_delim(samp_file[[i]],col_names = c("name", "filtered_reads", "value"),delim = " ") 
}

out <- bind_rows(samp) %>% dplyr::mutate(value=min(value)/value, group=mygroup) %>% dplyr::select(name, group, value)

write_tsv(out, outfile,col_names = F)
