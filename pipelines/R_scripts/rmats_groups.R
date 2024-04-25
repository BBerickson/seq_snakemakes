#!/usr/bin/env Rscript

library("tidyverse")

mydir <- snakemake@params[["project"]]
GRPS_UNIQ <- snakemake@params[["mygroup"]] %>% str_split_fixed(.,":",2) %>% c() %>% paste0("-",.,"$")
SAM_GRPS <- snakemake@params[["mysamps"]]


grps <- SAM_GRPS[str_detect(SAM_GRPS,GRPS_UNIQ[1])] %>% str_remove(.,GRPS_UNIQ[1]) %>% paste0(.,"*.bam$")
bams <- list.files(path = mydir, pattern = glob2rx(paste0("^", grps,collapse = "|"))) %>% paste0(mydir,"/",.)
# Combine pairs with commas
bams <- paste(bams, collapse = ",")
write_file(bams, snakemake@output[[1]])

grps <- SAM_GRPS[str_detect(SAM_GRPS,GRPS_UNIQ[2])] %>% str_remove(.,GRPS_UNIQ[2]) %>% paste0(.,"*.bam$")
bams <- list.files(path = mydir, pattern = glob2rx(paste0("^", grps,collapse = "|"))) %>% paste0(mydir,"/",.)
bams <- paste(bams, collapse = ",")
write_file(bams, snakemake@output[[2]])
