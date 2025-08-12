#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library("tidyverse")
library("jsonlite")

# arguments 
dedup_file <- args[1] #PROJ + "/stats/" + PROJ + "_clumpify.tsv" | PROJ + "/stats/" + PROJ + "_" + INDEX_MAP + "_UMI_dedup.tsv"
filter_file <- args[2] # PROJ + "/stats/" + PROJ + "_bbduk.tsv" | PROJ + "/stats/" + PROJ + "_cutadapt.tsv"
aligner_file <- args[3] # PROJ + "/stats/" + PROJ + "_aligned.tsv" 
FC_files <- strsplit(args[4], " ")[[1]] #expand( PROJ + "/stats/{sample}_summary_featureCounts.tsv",sample = SAMS_UNIQ)

index_map <- args[5] # INDEX_MAP,
sam_new <- args[6] # SAMPLES, 
mydir <- args[7] # PROJ
filter_type <- args[8] #  c("chrM|snRNA|snoRNA|rRNA|protein_coding")
outfile <- args[9]

sam_new <- fromJSON(gsub("'", '"', sam_new))
# file name and new name
sn <- enframe(sam_new, name = "new_name", value = "sample") %>%
  unnest(sample) %>%
  group_by(new_name) %>%
  mutate(new_name = ifelse(row_number() == 1, new_name, NA)) %>%
  ungroup()

# pre-alignment
if(str_detect(dedup_file,"clumpify")){
  dedup1 <- read_tsv(dedup_file,col_names = c("sample","type","value"),show_col_types = FALSE) %>% 
    spread(type,value) %>% 
    dplyr::mutate(Total_reads = Reads_In, 
                  Duplicates_removed = paste0(round(Duplicates_Found/Reads_In*100,digits = 2),"%")) %>% 
    dplyr::select(sample,Total_reads,Duplicates_removed)
  dedup2 <- read_tsv(filter_file,col_names = c("sample","type","value"),show_col_types = FALSE)%>% 
    separate(., "value", c("value","bbduk_reads_removed","empty"), "[\\(|\\)]") %>% 
    dplyr::select(sample,"bbduk_reads_removed")
  
  pre_alignment <- full_join(dedup1,dedup2,by="sample") 
  
} else {
  dedup1 <- read_tsv(filter_file,col_names = c("sample","type","value"),
                     show_col_types = FALSE) %>% 
    dplyr::mutate(type=str_replace_all(type," ","_")) %>% 
    spread(type,value) %>% 
    dplyr::mutate(Total_reads = as.double(Total_read_pairs_processed)) %>%
    dplyr::rename("passing_filters_Cutadapt"="Pairs_written_(passing_filters)") %>%
    dplyr::select(sample,Total_reads,passing_filters_Cutadapt)
  dedup2 <- read_tsv(dedup_file, col_names = c("sample", "type","value"), show_col_types = FALSE) %>% 
    dplyr::mutate(sample=str_remove(sample, paste0("_",index_map,"_UMI"))) %>%
    dplyr::mutate(type=str_replace_all(type," ","_")) 
  
  pre_alignment <- full_join(dedup1,dedup2,by="sample") %>% pivot_wider(names_from = type,values_from = value)
}

# alignment
aligner <- read_tsv(aligner_file,col_names = c("sample","type","alignment_rate"),show_col_types = FALSE) %>% 
  dplyr::filter(type == "overall_alignment_rate") %>% 
  dplyr::select(sample,alignment_rate)

# gather results
sn <- full_join(sn,aligner,by="sample") %>%
  full_join(.,pre_alignment,by="sample") %>% 
  arrange(new_name) %>% 
  mutate(new_name=if_else(is.na(new_name),sample,new_name))

# optional sub-sampling w/wo masking
mydirs <- paste0(mydir,"/bams_sub")
mydirm <- paste0(mydir,"/bams_mask")
mydirb <- paste0(mydir,"/bams")
mydirf <- paste0(mydir,"/stats")
if(!is_empty(list.files(path = mydirs, pattern = '.bam',recursive = T))){
  mask_count <- list.files(path = mydirm, pattern = '_mask_count',recursive = T)
  if(!is_empty(mask_count)){
    paste0(mydirm,"/",mask_count) -> mask_count
    
    gg <- NULL
    for(i in seq_along(mask_count)){
      gg <- read_delim(mask_count[i],delim = " ",
                       col_names = c("sample","index","value"),
                       show_col_types = FALSE) %>%
        bind_rows(gg)
      
    }
    sn <- gg %>% 
      separate(sample,into=c("sample","index"),sep="_(?!.*_)",extra="merge") %>% 
      dplyr::mutate(index = paste0("aligned_", index)) %>%
      pivot_wider(names_from = "index",values_from = c("value")) %>% 
      full_join(sn,.,by="sample")
    
    subsample_files <- list.files(path = mydirs, pattern = '_subsample.txt',recursive = T)
    if(!is_empty(subsample_files)){
      subsample <- NULL
      for(i in seq_along(subsample_files)){
        subsample <- read_delim(paste0(mydirs,"/",subsample_files[i]),
                                col_names = c("sample", "type","value"),
                                delim = " ",
                                show_col_types = FALSE) %>%
          dplyr::select(-type) %>% 
          bind_rows(subsample)
      }
      
      subsample_files <- list.files(path = mydirf, pattern = '_subsample_frac.tsv',recursive = T)
      
      if(!is_empty(subsample_files) & !is.null(subsample)){
        subsample <- read_tsv(paste0(mydirf,"/",subsample_files),
                              col_names = c("sample", "sub_group","read_fraction"),
                              show_col_types = FALSE) %>%
          full_join(subsample,.,by="sample") %>% 
          dplyr::mutate(read_fraction=paste0(value,"(", round(read_fraction*100,2),"%)")) %>%
          dplyr::select(-value) %>% 
          tidyr::separate(.,sample,c("sample","index"),"_(?=[^_]+$)") %>% 
          dplyr::mutate(index = paste0("subsampled_", index)) %>%
          pivot_wider(names_from = "index",values_from = c("read_fraction"))
        sn <- full_join(sn,subsample,by="sample") %>% arrange(sub_group,new_name)
      }
    }
    
  } else {
    subsample_files <- list.files(path = mydirs, pattern = '_subsample.txt',recursive = T)
    if(!is_empty(subsample_files)){
      subsample <- NULL
      for(i in seq_along(subsample_files)){
        subsample <- read_delim(paste0(mydirs,"/",subsample_files[i]),
                                col_names = c("sample", "type","value"),
                                delim = " ",
                                show_col_types = FALSE) %>%
          dplyr::select(-type) %>% 
          bind_rows(subsample)
      }
      subsample_files <- list.files(path = mydirf, pattern = '_subsample_frac.tsv',recursive = T)
      
      if(!is_empty(subsample_files) & !is.null(subsample)){
        subsample <- read_tsv(paste0(mydirf,"/",subsample_files),
                              col_names = c("sample", "sub_group","read_fraction"),
                              show_col_types = FALSE) %>%
          inner_join(subsample,.,by="sample") %>% 
          dplyr::mutate(read_fraction=paste0(value,"(", round(read_fraction*100,2),"%)")) %>% 
          separate(sample,into=c("sample","index"),sep="_(?!.*_)",extra="merge") %>% 
          dplyr::select(-value) %>% mutate(index=paste0("subsampled_", index)) %>% 
          pivot_wider(names_from = "index",values_from = "read_fraction")
        sn <- full_join(sn,subsample,by="sample") %>% arrange(sub_group,sample)
      }
    }
  }
}

# featureCounts
if(!is_empty(FC_files)){
  FC <- NULL
  for(i in seq_along(sn$sample)){
    FC <- read_delim(FC_files[str_detect(FC_files,sn$sample[i])],delim = " ",
                     col_names = c("type","value"),
                     show_col_types = FALSE) %>% 
      dplyr::mutate(sample=sn$sample[i]) %>% 
      dplyr::filter(!str_detect(type,filter_type)) %>%
      dplyr::mutate(value = if_else(str_detect(value, "spikin_"),1000000/value,value)) %>% 
      bind_rows(FC)
    
  }
  
  FC <- FC %>% distinct(type,sample,.keep_all =T) %>% spread(type,value) %>% 
    dplyr::select(sample, distinct(FC,type)$type) 
  
  sn <- full_join(sn,FC,by="sample") 
}

frag_files <- list.files(path = mydirf, pattern = "_fragment_results.tsv")

# fragment size
if(!is_empty(frag_files)){
  frag_files <- c(paste0(mydirf,"/",frag_files))
  frag <- frag_files %>%
    map(~read_tsv(.x, show_col_types = FALSE)) %>%
    reduce(full_join, by = "sample")
  
  sn <- full_join(sn,frag,by="sample") 
}

write_tsv(sn, outfile)
