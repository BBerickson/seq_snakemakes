#!/usr/bin/env Rscript

library("tidyverse")

clump <- read_tsv(snakemake@input[[1]],col_names = c("sample","type","value"),show_col_types = FALSE) 
bbduk <- read_tsv(snakemake@input[[2]],col_names = c("sample","type","value"),show_col_types = FALSE)
aligner <- read_tsv(snakemake@input[[3]],col_names = c("sample","type","alignment_rate"),show_col_types = FALSE)

mydir <- snakemake@params[["project"]]
cc <- clump %>% 
  spread(type,value) %>% 
  dplyr::mutate(Total_reads = Reads_In/2, 
                Duplicates_removed = paste0(round(Duplicates_Found/Reads_In*100,digits = 2),"%")) %>% 
  dplyr::select(sample,Total_reads,Duplicates_removed)

bb <- bbduk %>% 
  separate(., "value", c("value","bbduk_reads_removed","empty"), "[\\(|\\)]") %>% 
  dplyr::select(sample,"bbduk_reads_removed")

hh <- aligner %>% dplyr::select(sample,alignment_rate)

count_files <- list.files(path = mydir, pattern = '_count.txt',recursive = T)
count_files_names <- tibble(samp=count_files) %>% 
  dplyr::mutate(samp=str_remove(samp,"_count.txt")) %>% 
  unlist()

paste0(mydir,"/",count_files) -> count_files
gg <- NULL
for(i in seq_along(count_files)){
  if(str_detect(count_files[i],"_OR_count.txt")){
    gg <- read_delim(count_files[i],delim = " ",
                     col_names = c("sample", "type","value"),
                     show_col_types = FALSE) %>% 
      dplyr::filter(type != "NONE") %>% bind_rows(gg)
  } else {
    gg <- read_delim(count_files[i],delim = " ",
                     col_names = c("type","value"),
                     show_col_types = FALSE) %>% 
      dplyr::mutate(sample=count_files_names[i]) %>% 
      dplyr::filter(type != "NONE") %>% bind_rows(gg)
  }
  
}

gg <- gg %>% distinct(type,sample,.keep_all =T) %>% spread(type,value) %>% dplyr::select(sample, distinct(gg,type)$type)

full_join(cc,bb,by="sample") %>% 
  full_join(.,hh,by="sample") %>% 
  full_join(.,gg,by="sample") %>% 
  arrange(sample) %>% 
  write_tsv(., snakemake@output[[1]])

