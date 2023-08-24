#!/usr/bin/env Rscript

library("tidyverse")

# reads in bbduk stats file and mods to use reads removed #s
bbduk_kt <- read_lines(snakemake@input[[1]]) 
bbduk_fl <- read_lines(snakemake@input[[2]]) 
bbduk_fl <- bbduk_fl[str_detect(bbduk_fl,"Total Removed:")] %>% tibble(X1=.) %>% 
  separate(col = "X1",into=c("type","value","reads"),"\t") %>% select(value) %>% 
  separate(value,c("type","value","reads")," ",extra="drop") %>% select(-value) %>% 
  transmute(value=paste0(type,reads)) %>% str_replace(.,"\\(","\t") %>% str_remove(.,"\\)")

bbduk_kt[str_detect(bbduk_kt,"#Matched")] <- paste0("#Matched\t",bbduk_fl)

bbduk_out <- bbduk_kt
clumpfy <- read_lines(snakemake@input[[3]])
clumpfy <- clumpfy[str_detect(clumpfy,"Reads In:|Duplicates Found:")] %>% tibble(X1=.) %>% 
  separate(col = "X1",into=c("type","temp", "value"),"\\s",extra="merge") %>% 
  mutate(value=as.numeric(value)) %>% select(-temp) 
clumpfy2 <- clumpfy %>%
  # pivot_wider(names_from = type, values_from = value) %>%
  spread(type,value) %>% 
  transmute(value=paste0(Duplicates,"\t",Duplicates/Reads*100,"%"))
bbduk_kt[str_detect(bbduk_kt,"#Total")] <- paste0("#Total\t",clumpfy %>% filter(type=="Reads") %>% select(value))
bbduk_kt[str_detect(bbduk_kt,"#Matched")] <- paste0("#Matched\t",clumpfy2)
write_lines(bbduk_kt,snakemake@output[[2]])
# insure both files are made before next rule runs
write_lines(bbduk_out,snakemake@output[[1]])


