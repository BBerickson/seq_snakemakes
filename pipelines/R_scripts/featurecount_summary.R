#!/usr/bin/env Rscript

library("tidyverse")
# summaries featurecount output files made with GTF file sample and spikein
samp <- list()
col_names <- c("ID",	"Chr",	"Start", "End",	"Strand",	"Length", "gene_name",	"Geneid", "count")

for(i in seq_along(snakemake@input)) {
  index <- snakemake@params[["index"]][i]
  sub_reads <- read_tsv(snakemake@input[[i]],comment = "#",show_col_types = FALSE) 
  names(sub_reads) <- col_names
  sub_reads <- sub_reads %>% filter(Geneid %in% c("protein_coding","snRNA","snoRNA","rRNA")) %>% 
    mutate(Geneid = if_else(str_detect(Chr,"chrM|chrMT") & Geneid == "protein_coding", "chrM",Geneid))
  sub_reads <- sub_reads %>% 
    dplyr::select(Geneid, count) %>%
    group_by(Geneid) %>% summarise(count=sum(count,na.rm = T)) %>% 
    mutate(Geneid = paste0(Geneid,"_",index)) 
  myindex <- tibble(Geneid = index)
  samp[[i]] <- read_tsv(paste0(snakemake@input[[i]],".summary"),comment = "#",show_col_types = FALSE) %>% 
    dplyr::rename(count=names(.)[2]) %>% summarise(count=sum(count)) %>% bind_cols(myindex,.) %>% 
    bind_rows(.,sub_reads) %>% arrange(desc(count))
  
  if(file.exists(paste0(str_remove(snakemake@input[[i]],".tsv"),"_enrich.tsv"))){
    sub_reads <- read_tsv(paste0(str_remove(snakemake@input[[i]],".tsv"),"_enrich.tsv"),comment = "#",show_col_types = FALSE) 
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

write_delim(out, snakemake@output[[1]],col_names = F,delim = " ")
