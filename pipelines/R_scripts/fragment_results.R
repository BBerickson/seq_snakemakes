#!/usr/bin/env Rscript

library("tidyverse")

resl <- read_tsv(snakemake@input[[1]], show_col_types = FALSE) 
frag <- read_tsv(snakemake@input[[2]], show_col_types = FALSE)
if(length(snakemake@input)>2){
  frag2 <- read_tsv(snakemake@input[[3]], show_col_types = FALSE)
} else {
  frag2 <- NULL
}

samp <- paste(snakemake@params[["samp"]],c("Frag_Len_Mean", "Read_Len_Mean"),sep = "_")

names(frag)[1] <- "sample"
fragg <- frag %>% dplyr::select(sample,`Frag. Len. Mean`, `Read Len. Mean`) %>%
  dplyr::mutate(`Frag. Len. Mean` = round(`Frag. Len. Mean`,digits = 5),
         `Read Len. Mean` = round(`Read Len. Mean`,digits = 5)) %>% 
  dplyr::rename(!!samp[1] := `Frag. Len. Mean`, !!samp[2] := `Read Len. Mean`)

if(!is.null(frag2)){
  spik <-paste(snakemake@params[["spik"]],c("Frag_Len_Mean", "Read_Len_Mean"),sep = "_")
  names(frag2)[1] <- "sample"
  fragg <- frag2 %>% dplyr::select(sample,`Frag. Len. Mean`, `Read Len. Mean`) %>%
    dplyr::mutate(`Frag. Len. Mean` = round(`Frag. Len. Mean`,digits = 5),
           `Read Len. Mean` = round(`Read Len. Mean`,digits = 5)) %>%
    dplyr::rename(!!spik[1] := `Frag. Len. Mean`, !!spik[2] := `Read Len. Mean`) %>% 
    full_join(fragg,.,by="sample")
}

full_join(resl,fragg,by="sample") %>%
  arrange(sample) %>% 
  write_tsv(., snakemake@input[[1]])

write_tsv(fragg, snakemake@output[[1]])
