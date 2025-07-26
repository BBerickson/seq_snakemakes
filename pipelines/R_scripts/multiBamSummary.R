#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("tidyverse")

scalefactor <- read_tsv(args[1],show_col_types = FALSE) 
counts <- args[2]
counts <- strsplit(counts," ")[[1]]
index <- args[3]

for(i in scalefactor$sample){
  cc <- counts[str_detect(counts,paste0(i,"_summary_featureCounts.tsv$"))]
  ccc <- read_delim(cc,col_names = c("sample","scalingFactor"),delim = " ",show_col_types = FALSE)
  
  ccc %>% bind_rows(.,scalefactor %>% filter(sample == i) %>% mutate(sample=paste0("spikin_",index[1]),scalingFactor=1000000/scalingFactor)) %>%
    distinct() %>% 
    write_delim(., cc,col_names = F,delim = " ")
}


