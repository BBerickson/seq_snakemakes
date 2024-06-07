#!/usr/bin/env Rscript

library("tidyverse")

samp <- paste(snakemake@params[["samp"]],c("Frag_Len_Mean", "Read_Len_Mean"),sep = "_")
resl <- read_tsv(snakemake@input[[1]], show_col_types = FALSE) 

mydir <- snakemake@params[["project"]]
mydirf <- paste0(mydir,"/stats")

my_files <- list.files(path = mydirf, pattern = '_fragment.txt',recursive = T) 
my_files1 <- paste0(mydirf,"/",my_files[str_detect(my_files,pattern = snakemake@params[["samp"]])])

if(length(my_files) > length(my_files1)){
  my_files2 <- paste0(mydirf,"/",my_files[str_detect(my_files,pattern = snakemake@params[["spik"]])])
} else {
  my_files2 <- NULL
}


frag <- NULL
for(i in my_files1){
  frag <- read_tsv(i, show_col_types = FALSE) %>% 
    dplyr::rename(sample=`...1`) %>% 
    dplyr::select(sample,`Frag. Len. Mean`, `Read Len. Mean`) %>%
    dplyr::mutate(`Frag. Len. Mean` = round(`Frag. Len. Mean`,digits = 5),
                  `Read Len. Mean` = round(`Read Len. Mean`,digits = 5)) %>% 
    dplyr::rename(!!samp[1] := `Frag. Len. Mean`, !!samp[2] := `Read Len. Mean`) %>% 
    bind_rows(frag,.)
}

frag <- frag %>% separate(sample,into=c("dir","sample"),sep="/bams/") %>% 
  separate(sample,into=c("sample","file"),sep=paste0("_", snakemake@params[["samp"]])) %>% 
  dplyr::select(-dir,-file)


if(!is.null(my_files2)){
  frag2 <- NULL
  spik <-paste(snakemake@params[["spik"]],c("Frag_Len_Mean", "Read_Len_Mean"),sep = "_")
  for(i in my_files2){
    frag2 <- read_tsv(i, show_col_types = FALSE) %>% 
      dplyr::rename(sample=`...1`) %>% 
      dplyr::select(sample,`Frag. Len. Mean`, `Read Len. Mean`) %>%
      dplyr::mutate(`Frag. Len. Mean` = round(`Frag. Len. Mean`,digits = 5),
                    `Read Len. Mean` = round(`Read Len. Mean`,digits = 5)) %>% 
      dplyr::rename(!!spik[1] := `Frag. Len. Mean`, !!spik[2] := `Read Len. Mean`) %>% 
      bind_rows(frag2,.)
  }
  frag2 <- frag2 %>% separate(sample,into=c("dir","sample"),sep="/bams/") %>% 
    separate(sample,into=c("sample","file"),sep=paste0("_", snakemake@params[["spik"]])) %>% 
    dplyr::select(-dir,-file)
  
  frag <- full_join(frag2, frag,.,by="sample")
} 

full_join(resl,frag,by="sample") %>%
  arrange(sample) %>% 
  write_tsv(., snakemake@input[[1]])

write_tsv(frag, snakemake@output[[1]])
