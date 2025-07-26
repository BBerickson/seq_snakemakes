#!/usr/bin/env Rscript

library("tidyverse")

clump <- read_tsv(snakemake@input[[1]],col_names = c("sample","type","value"),show_col_types = FALSE) 
bbduk <- read_tsv(snakemake@input[[2]],col_names = c("sample","type","value"),show_col_types = FALSE)
aligner <- read_tsv(snakemake@input[[3]],col_names = c("sample","type","alignment_rate"),show_col_types = FALSE)

mydir <- snakemake@params[["project"]]
mydirc <- paste0(mydir,"/stats")
mydirf <- paste0(mydir,"/stats")
sam_new <- snakemake@params[["sam_new"]]

sn <- tibble(
  new_name = names(sam_new),
  sample = unlist(sam_new)
)

cc <- clump %>% 
  spread(type,value) %>% 
  dplyr::mutate(Total_reads = Reads_In, 
                Duplicates_removed = paste0(round(Duplicates_Found/Reads_In*100,digits = 2),"%")) %>% 
  dplyr::select(sample,Total_reads,Duplicates_removed)

bb <- bbduk %>% 
  separate(., "value", c("value","bbduk_reads_removed","empty"), "[\\(|\\)]") %>% 
  dplyr::select(sample,"bbduk_reads_removed")

hh <- aligner %>% 
  dplyr::filter(type == "overall_alignment_rate") %>% 
  dplyr::select(sample,alignment_rate)

count_files <- list.files(path = mydirc, pattern = '_summary_featureCounts.tsv',recursive = T)
if(!is_empty(count_files)){
  count_files_names <- tibble(samp=count_files) %>% 
    dplyr::mutate(samp=str_remove(samp,"_summary_featureCounts.tsv")) %>%
    unlist()
  
  paste0(mydirc,"/",count_files) -> count_files
  filter_type <- c("NONE")
  
  gg <- NULL
  for(i in seq_along(count_files)){
    gg <- read_delim(count_files[i],delim = " ",
                     col_names = c("type","value"),
                     show_col_types = FALSE) %>% 
      dplyr::mutate(sample=count_files_names[i]) %>% 
      dplyr::filter(!str_detect(type,filter_type)) %>%
      dplyr::mutate(value = if_else(str_detect(value, "spikin_"),1000000/value,value)) %>% 
      bind_rows(gg)
    
  }
  
  
  hh <- gg %>% distinct(type,sample,.keep_all =T) %>% spread(type,value) %>% 
    dplyr::select(sample, distinct(gg,type)$type) %>% 
    full_join(hh,.,by="sample")
}

frag_files <- list.files(path = mydirf, pattern = '_fragment_results.tsv',recursive = T)
if(!is_empty(frag_files)){
  frag_files <- paste0(mydirf,"/",frag_files)
  gg <- NULL
  for(i in seq_along(frag_files)){
    gg <- read_tsv(frag_files[i],
                   show_col_types = FALSE) %>% 
      mutate(sample=str_remove(sample,"_spike")) %>% 
      pivot_longer(cols = !sample,names_to = "type",values_to = "length") %>% 
      bind_rows(gg)
    
  }
  
  hh <- gg %>% pivot_wider(names_from = type,values_from = length) %>% 
    full_join(hh,.,by="sample")
}


full_join(sn,cc,by="sample") %>% 
  full_join(.,bb,by="sample") %>% 
  full_join(.,hh,by="sample") %>% 
  arrange(new_name) %>% 
  write_tsv(., snakemake@output[[1]])

