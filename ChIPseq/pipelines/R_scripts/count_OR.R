#!/usr/bin/env Rscript

library("tidyverse")

type <- snakemake@params[["type"]]
samp <- read_delim(snakemake@input[[1]],col_names = c("name","value"),delim = " ") 
if(length(snakemake@input)>1){
  input <- read_delim(snakemake@input[[2]],col_names = c("name","value"),delim = " ")
} else {
  input <- NULL
}

if(str_detect(type,"OR|OR_enrich") & !is.null(input)){
  sam_en <- snakemake@params[["sam_en"]]
  spk_en <- snakemake@params[["spk_en"]]
  in_sam <- snakemake@params[["in_sam"]]
  in_spk <- snakemake@params[["in_spk"]]
  InputSpike <- filter(input,str_detect(name,paste0("^",in_spk,"$"))) %>% dplyr::rename(value1=value) %>% mutate(name="test")
  ChIP_en <- filter(samp,str_detect(name,paste0("^",sam_en,"$"))) %>% dplyr::rename(value2_en=value) %>% mutate(name="test")
  ChIP <- filter(samp,str_detect(name,paste0("^",in_sam,"$"))) %>% dplyr::rename(value2=value) %>% mutate(name="test")
  Input <- filter(input,str_detect(name,paste0("^",in_sam,"$"))) %>% dplyr::rename(value3=value) %>% mutate(name="test")
  ChIPSpike_en <- filter(samp,str_detect(name,paste0("^",spk_en,"$"))) %>% dplyr::rename(value4_en=value) %>% mutate(name="test")
  ChIPSpike <- filter(samp,str_detect(name,paste0("^",in_spk,"$"))) %>% dplyr::rename(value4=value) %>% mutate(name="test")
  
  # OR = (InputSpike*ChIP)/(Input*ChIPSpike)
  out <- full_join(ChIPSpike,ChIP,by="name") %>% 
    full_join(Input,by="name") %>% 
    full_join(InputSpike,by="name") %>% 
    summarise(value=round((value1*value2)/(value3*value4),digits = 5)) %>% 
    mutate(name= "test", type= "OR")
  out <- full_join(ChIPSpike_en,ChIP_en,by="name") %>% 
    full_join(Input,by="name") %>% 
    full_join(InputSpike,by="name") %>% 
    summarise(value=round((value1*value2_en)/(value3*value4_en),digits = 5)) %>% 
    mutate(name= "test", type= "OR_enrich") %>% 
    bind_rows(out,.)
  out1 <- out %>% full_join(ChIP,.,by="name") %>% 
    mutate(value=round(10000000/as.integer(value2/value),digits = 5), type=paste0(type,"_",in_sam)) %>% 
    dplyr::select(-value2) 
  out <- out %>% full_join(ChIPSpike,.,by="name") %>% 
    mutate(value=round(10000000/as.integer(value4/value),digits = 5), type=paste0(type,"_",in_spk)) %>% 
    dplyr::select(-value4) %>% 
    bind_rows(out,out1,.)
} else {
  ChIP <- filter(samp,str_detect(name,paste0("^",type,"$"))) %>% mutate(name="test","type"=type)
  out <- tibble(name="test", value = 1e6, type=type) %>% 
    bind_rows(.,ChIP)
}

out <- out %>% 
  dplyr::mutate(mysamp=snakemake@input[[1]]) %>% 
  separate(.,mysamp,c("temp","mysamp"),"/bams/") %>% 
  mutate(mysamp=str_remove(mysamp,"_count.txt")) %>% 
  select(mysamp,type,value)
write_delim(out, snakemake@output[[1]],col_names = F,delim = " ")
