#!/usr/bin/env Rscript

library("tidyverse")

clump <- read_tsv(snakemake@input[[1]],col_names = c("sample","type","value"),show_col_types = FALSE) 
bbduk <- read_tsv(snakemake@input[[2]],col_names = c("sample","type","value"),show_col_types = FALSE)
aligner <- read_tsv(snakemake@input[[3]],col_names = c("sample","type","alignment_rate"),show_col_types = FALSE)

cc <- clump %>% 
  spread(type,value) %>% 
  dplyr::mutate(Total_reads = Reads_In/2, 
                Duplicates_removed = paste0(round(Duplicates_Found/Reads_In*100,digits = 2),"%")) %>% 
  dplyr::select(sample,Total_reads,Duplicates_removed)

bb <- bbduk %>% 
  separate(., "value", c("value","bbduk_reads_removed","empty"), "[\\(|\\)]") %>% 
  dplyr::select(sample,"bbduk_reads_removed")

hh <- aligner %>%  dplyr::filter(type == "overall alignment rate") %>% 
  mutate(type=str_replace_all(type," ","_")) %>% dplyr::select(sample,alignment_rate)

mydir <- snakemake@params[["project"]]


subsample_files <- list.files(path = mydir, pattern = '_subsample.tsv',recursive = T)
subsample <- NULL
for(i in seq_along(subsample_files)){
  subsample <- read_delim(paste0(mydir,"/",subsample_files[i]),
                          col_names = c("sample", "type","value"),
                          delim = " ",
                          show_col_types = FALSE) %>% 
    mutate(type=str_replace_all(type," ","_")) %>%
    separate(sample,into=c("sample","index"),sep="_(?!.*_)",extra="merge") %>% 
    bind_rows(subsample)
}

subsample <- subsample %>% 
  spread(type,value) 

subsample_files <- list.files(path = mydir, pattern = '_subsample_frac.tsv',recursive = T)
subsample <- read_tsv(paste0(mydir,"/",subsample_files[1]),
                      col_names = c("sample", "group","read_fraction"),
                      show_col_types = FALSE) %>%
  dplyr::mutate(read_fraction=round(read_fraction,3)) %>% 
  separate(sample,into=c("sample","index"),sep="_(?!.*_)",extra="merge") %>% 
  full_join(subsample,.,by=c("sample","index")) %>% 
  pivot_wider(names_from = "index",values_from = c("Filtered_reads","read_fraction","Sampled_reads"))%>% 
  dplyr::select(sample,starts_with("Filtered_reads"),group,starts_with("read_fraction"),starts_with("Sampled_reads")) 

full_join(cc,bb,by="sample") %>% 
  full_join(.,hh,by="sample") %>% 
  full_join(.,subsample,by="sample") %>% 
  arrange(sample) %>% 
  write_tsv(., snakemake@output[[1]])

