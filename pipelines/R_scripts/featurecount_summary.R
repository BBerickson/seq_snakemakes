#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library("tidyverse")
# summaries featurecount output files made with GTF file sample and spikein

featurecount_files <- strsplit(args[1], " ")[[1]]
indexs <- strsplit(args[2], " ")[[1]]
outfile <- args[3]

samp <- list()
col_names <- c("ID",	"Chr",	"Start", "End",	"Strand",	"Length", "gene_name",	"Geneid", "count")

for(i in seq_along(featurecount_files)) {
  index <- indexs[str_detect(basename(featurecount_files[[i]]),indexs)][1]
  sub_reads <- read_tsv(featurecount_files[[i]],comment = "#",show_col_types = FALSE) 
  names(sub_reads) <- col_names
  sub_reads_filter <- sub_reads %>% filter(Geneid %in% c("protein_coding","snRNA","snoRNA","rRNA"))
  if(nrow(sub_reads_filter)==0){
    sub_reads_filter <- sub_reads
  } 
    
  sub_reads_filter <- sub_reads_filter %>% 
    mutate(Geneid = if_else(str_detect(Chr,"chrM|chrMT") & Geneid == "protein_coding", "chrM",Geneid))
  sub_reads_filter <- sub_reads_filter %>% 
    dplyr::select(Geneid, count) %>%
    mutate(Geneid=if_else(is.na(Geneid),"Genes",Geneid)) %>% 
    group_by(Geneid) %>% summarise(count=sum(count,na.rm = T)) %>% 
    mutate(Geneid = paste0(Geneid,"_",index)) 
  myindex <- tibble(Geneid = index)
  samp[[i]] <- read_tsv(paste0(featurecount_files[[i]],".summary"),comment = "#",show_col_types = FALSE) %>% 
    dplyr::rename(count=names(.)[2]) %>% summarise(count=sum(count)) %>% bind_cols(myindex,.) %>% 
    bind_rows(.,sub_reads_filter) %>% arrange(desc(count))
  
  if(file.exists(paste0(str_remove(featurecount_files[[i]],".tsv"),"_enrich.tsv"))){
    sub_reads <- read_tsv(paste0(str_remove(featurecount_files[[i]],".tsv"),"_enrich.tsv"),comment = "#",show_col_types = FALSE) 
    names(sub_reads) <- col_names
    sub_reads <- sub_reads %>% 
      mutate(Geneid = if_else(str_detect(Chr,"chrM|chrMT"), Chr,Geneid)) %>% 
      filter(Geneid == "protein_coding") %>% 
      dplyr::select(Geneid, count) %>%
      group_by(Geneid) %>% summarise(count=sum(count,na.rm = T)) %>% 
      mutate(Geneid = paste0("enrich_", index))
    samp[[i]] <- bind_rows(samp[[i]],sub_reads)  %>% arrange(desc(count))
  }
  
}

out <- bind_rows(samp)

write_delim(out, outfile,col_names = F,delim = " ")
