#!/usr/bin/env Rscript

library("tidyverse")

samp <- paste(snakemake@params[["samp"]],c("Frag_Len_Mean", "Read_Len_Mean"),sep = "_")

frag <- NULL
for(i in seq_along(snakemake@input)){
  frag_temp <- read_tsv(snakemake@input[[i]], show_col_types = FALSE) %>% 
    dplyr::rename(sample=`...1`) %>% 
    dplyr::select(sample,`Frag. Len. Mean`, `Read Len. Mean`) %>%
    dplyr::mutate(`Frag. Len. Mean` = round(`Frag. Len. Mean`,digits = 5),
                  `Read Len. Mean` = round(`Read Len. Mean`,digits = 5)) %>% 
    dplyr::rename(!!samp[1] := `Frag. Len. Mean`, !!samp[2] := `Read Len. Mean`) 
  
  if(str_detect(frag_temp$sample,"/bams/")){
    frag_temp <- frag_temp %>% separate(sample,into=c("dir","sample"),sep="/bams/") %>%
      separate(sample,into=c("sample","file"),sep=paste0("_", snakemake@params[["samp"]])) %>%
      dplyr::select(-dir,-file)
  }
  frag <- bind_rows(frag,frag_temp)
}


write_tsv(frag, snakemake@output[[1]])
