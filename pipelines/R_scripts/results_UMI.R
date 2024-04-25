#!/usr/bin/env Rscript

library("tidyverse")

cutad <- read_tsv(snakemake@input[[1]],col_names = c("sample","type","value"),
                  show_col_types = FALSE) 
aligner <- read_tsv(snakemake@input[[2]],col_names = c("sample","type","alignment_rate"),
                    show_col_types = FALSE) %>% 
  dplyr::filter(type == "overall_alignment_rate") %>% 
  dplyr::select(sample,alignment_rate)

mydir <- snakemake@params[["project"]]
mydirc <- paste0(mydir,"/counts")
mydirf <- paste0(mydir,"/stats")
mydirs <- paste0(mydir,"/bams_sub")
mydirb <- paste0(mydir,"/bams")
index_map <- snakemake@params[["index_map"]]

cutad <- cutad %>% 
  dplyr::mutate(type=str_replace_all(type," ","_")) %>% 
  spread(type,value) %>% 
  dplyr::mutate(Total_reads = as.double(Total_read_pairs_processed)) %>%
  dplyr::rename("passing_filters_Cutadapt"="Pairs_written_(passing_filters)") %>%
  dplyr::select(sample,Total_reads,passing_filters_Cutadapt)

dedup_files <- list.files(path = mydirf, pattern = '_dedup.tsv',recursive = T)
dedup <- read_tsv(paste0(mydirf,"/",dedup_files),
                    col_names = c("sample", "type","value"),
                    show_col_types = FALSE) %>% 
  dplyr::mutate(sample=str_remove(sample, paste0("_",index_map,"_UMI")),
                  index = index_map) %>%
  dplyr::mutate(type=str_replace_all(type," ","_")) 

dedup <- dedup %>% dplyr::filter(type == "Input_Reads" | type == "Number_of_reads_out" ) %>%
  spread(type,value) %>% dplyr::mutate(Unique_UMI=paste0(Number_of_reads_out,"(",paste0(round(Number_of_reads_out/Input_Reads*100,1),"%"),")")) %>% 
  full_join(aligner,.,by="sample") %>% 
  dplyr::mutate(Aligned_reads=paste0(Input_Reads,"(",alignment_rate,")")) %>% 
  dplyr::select(-Number_of_reads_out,-alignment_rate,-Input_Reads) %>% 
  pivot_wider(names_from = "index",values_from = c("Aligned_reads","Unique_UMI"))

hh <- full_join(cutad,dedup,by="sample") %>% 
  arrange(sample)

# tests if subsample per each genome or per mixed genome
count_files <- list.files(path = mydirb, pattern = '_filtred_count',recursive = T)
if(!is_empty(count_files)){
  
  paste0(mydirb,"/",count_files) -> count_files
  
  gg <- NULL
  for(i in seq_along(count_files)){
    gg <- read_delim(count_files[i],delim = " ",
                     col_names = c("sample","index","value"),
                     show_col_types = FALSE) %>%
      bind_rows(gg)
    
  }
  hh <- gg %>% 
    separate(sample,into=c("sample","index"),sep="_(?!.*_)",extra="merge") %>% 
    dplyr::mutate(index = paste0("aligned_", index)) %>%
    pivot_wider(names_from = "index",values_from = c("value")) %>% 
    full_join(hh,.,by="sample")
  
  subsample_files <- list.files(path = mydirs, pattern = '_subsample.txt',recursive = T)
  if(!is_empty(subsample_files)){
    subsample <- NULL
    for(i in seq_along(subsample_files)){
      subsample <- read_delim(paste0(mydirs,"/",subsample_files[i]),
                              col_names = c("sample", "type","value"),
                              delim = " ",
                              show_col_types = FALSE) %>%
        dplyr::select(-type) %>% 
        bind_rows(subsample)
    }
    
    subsample_files <- list.files(path = mydirf, pattern = '_subsample_frac.tsv',recursive = T)
    
    if(!is_empty(subsample_files) & !is.null(subsample)){
      subsample <- read_tsv(paste0(mydirf,"/",subsample_files),
                            col_names = c("sample", "sub_group","read_fraction"),
                            show_col_types = FALSE) %>%
        full_join(subsample,.,by="sample") %>% 
        dplyr::mutate(read_fraction=paste0(value,"(", round(read_fraction*100,2),"%)")) %>%
        dplyr::select(-value) %>% 
        tidyr::separate(.,sample,c("sample","index"),"_(?=[^_]+$)") %>% 
        dplyr::mutate(index = paste0("subsampled_", index)) %>%
        pivot_wider(names_from = "index",values_from = c("read_fraction"))
      hh <- full_join(hh,subsample,by="sample") %>% arrange(sub_group,sample)
    }
  }
  subsample_files <- list.files(path = mydirf, pattern = '_subsample.tsv',recursive = T)
  if(!is_empty(subsample_files)){
    subsample <- NULL
    for(i in seq_along(subsample_files)){
      subsample <- read_delim(paste0(mydirf,"/",subsample_files[i]),
                              col_names = c("sample", "type","value"),
                              delim = " ",
                              show_col_types = FALSE) %>%
        mutate(type=str_replace_all(type," ","_")) %>%
        separate(sample,into=c("sample","index"),sep="_(?!.*_)",extra="merge") %>%
        bind_rows(subsample)
    }
    
    subsample <- subsample %>%
      dplyr::mutate(type = paste(type,index,sep="_")) %>% 
      dplyr::select(-index) %>% 
      spread(type,value)
    hh <- full_join(hh,subsample,by="sample")
  }
  
} else {
  subsample_files <- list.files(path = mydirc, pattern = '_subsample.txt',recursive = T)
  if(!is_empty(subsample_files)){
    subsample <- NULL
    for(i in seq_along(subsample_files)){
      subsample <- read_delim(paste0(mydirc,"/",subsample_files[i]),
                              col_names = c("sample", "type","value"),
                              delim = " ",
                              show_col_types = FALSE) %>%
        dplyr::select(-type) %>% 
        bind_rows(subsample)
    }
    subsample_files <- list.files(path = mydirf, pattern = '_subsample_frac.tsv',recursive = T)
    
    if(!is_empty(subsample_files) & !is.null(subsample)){
      subsample <- read_tsv(paste0(mydirf,"/",subsample_files),
                            col_names = c("sample", "sub_group","read_fraction"),
                            show_col_types = FALSE) %>%
        full_join(subsample,.,by="sample") %>% 
        dplyr::mutate(read_fraction=paste0(value,"(", round(read_fraction*100,2),"%)")) %>% 
        dplyr::mutate(sample=str_remove(sample, paste0("_",index_map))) %>%
        dplyr::rename(!!paste0("subsampled_", index_map) := read_fraction) %>% 
        dplyr::select(-value) 
      hh <- full_join(hh,subsample,by="sample") %>% arrange(sub_group,sample)
    }
  }
  
  count_files <- list.files(path = mydirc, pattern = '_count',recursive = T)
  if(!is_empty(count_files)){
    count_files_names <- tibble(samp=count_files) %>% 
      dplyr::mutate(samp=str_remove(samp,"_count.txt")) %>%
      dplyr::mutate(samp=str_remove(samp,"_spikin")) %>% 
      unlist()
    
    paste0(mydirc,"/",count_files) -> count_files
    
    gg <- NULL
    for(i in seq_along(count_files)){
      gg <- read_delim(count_files[i],delim = " ",
                       col_names = c("type","value"),
                       show_col_types = FALSE) %>% 
        dplyr::filter(!str_detect(type,"enrich_|chrM_")) %>% 
        dplyr::mutate(sample=count_files_names[i]) %>% 
        bind_rows(gg)
    }
    
    hh <- gg %>% distinct(type,sample,.keep_all =T) %>% spread(type,value) %>% 
      dplyr::select(sample, distinct(gg,type)$type) %>% 
      full_join(hh,.,by="sample")
  }
}



hh %>% 
  write_tsv(., snakemake@output[[1]])

