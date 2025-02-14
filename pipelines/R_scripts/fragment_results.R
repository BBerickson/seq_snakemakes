#!/usr/bin/env Rscript

library("tidyverse")

num_cols <- count_fields(snakemake@input[[1]],
                         n_max = 1,
             tokenizer = tokenizer_tsv())
frag <- NULL
if(num_cols > 2){
  samp <- paste(snakemake@params[["samp"]],c("Frag_Len_Mean", "Read_Len_Mean"),sep = "_")
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
} else {
  samp <- paste(snakemake@params[["samp"]],c("Read_Len_Mean", "Read_Len_Median"),sep = "_")
  for(i in seq_along(snakemake@input)){
    samp_name <- snakemake@input[[i]] %>% str_replace(., ".*/", "") %>% 
      str_remove(.,paste0("_", snakemake@params[["samp"]],"_fragment.txt"))
    
    frag_temp <- read_tsv(snakemake@input[[i]], show_col_types = FALSE, col_names = F, skip = 4, n_max = 2) %>% 
      mutate(sample=samp_name) %>% pivot_wider(values_from =X2,names_from = X1) %>% 
      dplyr::rename(!!samp[1] := `#Avg:`, !!samp[2] := `#Median:`)
    frag <- bind_rows(frag,frag_temp)
  }
}

write_tsv(frag, snakemake@output[[1]])
