#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library("tidyverse")

samp_file <- strsplit(args[1], " ")[[1]]
type <- args[2]
indexs <- strsplit(args[3], " ")[[1]]
outfile <- args[4]


samp <- read_delim(samp_file[[1]],col_names = c("name","value"),delim = " ",show_col_types = FALSE) 
if(length(samp_file)>1){
  input <- read_delim(samp_file[[2]],col_names = c("name","value"),delim = " ",show_col_types = FALSE)
} else {
  input <- NULL
}

if(!is.null(input) & length(indexs) > 1){
  in_sam <- indexs[[1]]
  in_spk <- indexs[[2]]
  sam_en <- paste0("enrich_",in_sam)
  spk_en <- paste0("enrich_",in_spk)
  InputSpike <- filter(input,str_detect(name,paste0("^",in_spk,"$"))) %>% dplyr::rename(value1=value) %>% mutate(name="test")
  ChIP_en <- filter(samp,str_detect(name,paste0("^",sam_en,"$"))) %>% dplyr::rename(value2_en=value) %>% mutate(name="test")
  ChIP <- filter(samp,str_detect(name,paste0("^",in_sam,"$"))) %>% dplyr::rename(value2=value) %>% mutate(name="test")
  Input <- filter(input,str_detect(name,paste0("^",in_sam,"$"))) %>% dplyr::rename(value3=value) %>% mutate(name="test")
  ChIPSpike_en <- filter(samp,str_detect(name,paste0("^",spk_en,"$"))) %>% dplyr::rename(value4_en=value) %>% mutate(name="test")
  ChIPSpike <- filter(samp,str_detect(name,paste0("^",in_spk,"$"))) %>% dplyr::rename(value4=value) %>% mutate(name="test")
  
  # OR = (InputSpike*ChIP)/(Input*ChIPSpike)
  out <- full_join(ChIPSpike_en,ChIP_en,by="name") %>% 
    full_join(Input,by="name") %>% 
    full_join(InputSpike,by="name") %>% 
    summarise(value=round((value1*value2_en)/(value3*value4_en),digits = 5)) %>% 
    mutate(name= "test", type= "OR_enrich") 
  out1 <- out %>% full_join(ChIP,.,by="name") %>% 
    mutate(value=as.integer(value*value2)) %>% 
    dplyr::select(-value2) %>% 
    dplyr::mutate(type=paste0(type,"_",in_sam)) 
  out <- out %>% full_join(ChIPSpike,.,by="name") %>% 
    mutate(value=as.integer(value*value4)) %>% 
    dplyr::select(-value4) %>% 
    dplyr::mutate(type=paste0(type,"_",in_spk)) %>%
    bind_rows(out,out1,.)
  out <- out %>% 
    dplyr::mutate(mysamp=samp_file[[1]]) %>% 
    separate(.,mysamp,c("temp","mysamp"),"/.*/") %>% 
    mutate(mysamp=str_remove(mysamp,"_summary_featureCounts.tsv")) %>%  
    select(mysamp,type,value)
  
  write_delim(out %>% select(-mysamp), samp_file[[1]],col_names = F,delim = " ", append = TRUE)
  write_delim(out, outfile,col_names = F,delim = " ")
} else {
  out <- tibble(mysamp=samp_file[[1]], type="OR", value = NA) %>% 
    separate(.,mysamp,c("temp","mysamp"),"/.*/") %>% 
    mutate(mysamp=str_remove(mysamp,"_summary_featureCounts.tsv")) %>%  
    select(mysamp,type,value)
  write_delim(out, outfile,col_names = F,delim = " ")
}


