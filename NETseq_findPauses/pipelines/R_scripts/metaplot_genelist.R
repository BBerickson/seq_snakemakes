#!/usr/bin/env Rscript

# make master in common gene lists
library(tidyverse)

gene_coords <- unlist(snakemake@input)
groups <- snakemake@params[["groups"]]
# extract the first file in each group, should correspond to the control file
first_items <- sapply(groups, function(x) x[1]) %>% paste0(.,collapse = "|")
gene_coords <- gene_coords[str_detect(gene_coords,first_items)]

beds <- list()
for(i in gene_coords){
  beds[[i]] <- read_delim(i,delim = "\t", col_names =c(
    "chrom", "start", 
    "end",   "name", 
    "score", "strand", "counts"
  ), show_col_types = FALSE) %>% mutate(score=i) %>% 
    distinct(name,.keep_all = T)
}

beds <- bind_rows(beds) %>% group_by(name) %>%
  filter(n_distinct(score) == length(gene_coords)) %>%
  ungroup() %>% select(-score) %>% mutate(gene_len = (end - start) / 1000) %>%
  group_by(chrom, start, end, name, strand, gene_len) %>%
  summarize(counts = mean(counts), .groups = "drop")

#### save file
write_tsv(beds, snakemake@output[[1]])
