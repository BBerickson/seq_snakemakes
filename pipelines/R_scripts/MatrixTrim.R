#!/usr/bin/env Rscript

# edit deeptools computeMatrix files
# marks overlapping genes based on strand, 
# filters out overlapping genes on same strand, filters out based on less than size
# outputs a zipped "matirx" file

# -options max_gene_length, gene_seperation, 

library("tidyverse")
library("valr")

my_args <- snakemake@params[["args"]] %>% 
  as_tibble() %>% 
  separate(.,value,c("size_select","dist_select","stranded")," ")

# meta data
meta <- read_lines(snakemake@input[[1]],n_max = 1)

num_bins <-
  count_fields(snakemake@input[[1]],
               n_max = 1,
               skip = 1,
               tokenizer = tokenizer_tsv())

tablefile <- suppressMessages(read_tsv(
  snakemake@input[[1]],
  comment = "@",
  col_names = c("chrom", "start", "end","gene", "value", "sign", 1:(num_bins - 6)))) 
num_genes <- n_distinct(tablefile$gene)

if(as.numeric(my_args$dist_select) > 0 & as.numeric(my_args$size_select) > 0){
  
  tablefile <- tablefile %>%  
    dplyr::mutate(value2 = end-start) 
  if(my_args$stranded == "stranded"){
    tablefile <- tablefile %>% 
      group_by(sign) %>% 
      bed_cluster(max_dist = as.numeric(my_args$dist_select)) %>% 
      ungroup() 
  } else {
    tablefile <- tablefile %>%
      bed_cluster(max_dist = as.numeric(my_args$dist_select))
  }
  
  tablefile <- tablefile %>% 
    arrange(.id) %>%
    dplyr::filter(!(duplicated(.id) | duplicated(.id, fromLast=TRUE))) %>%
    dplyr::filter(value2 >= as.numeric(my_args$size_select)) %>% 
    dplyr::select(-value2, -.id)
}

meta <- meta %>% str_replace(paste0("\\[0,",num_genes,"\\]"),paste0("\\[0,",nrow(tablefile),"\\]"))

write_lines(meta,snakemake@output[[1]])
write_tsv(tablefile, snakemake@output[[1]], col_names = F,append = T)

