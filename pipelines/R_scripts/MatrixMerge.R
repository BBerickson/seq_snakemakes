#!/usr/bin/env Rscript

# merge deeptools computeMatrix files 

library("tidyverse")

my_names <- snakemake@params[["samplab"]] 
groulab <- snakemake@params[["groulab"]] 
my_names <- my_names[[groulab]]

my_files <- unlist(snakemake@input)
sense_files <- my_files[!str_detect(my_files,"_anti_matrix.gz")]
asense_files <- my_files[str_detect(my_files,"_anti_matrix.gz")]

if(!length(sense_files)>0) {
  sense_files <- asense_files
  asense_files <- NULL
  my_names <- str_c(my_names,"_antisense")
} 

num_bins <-
  count_fields(sense_files[1],
               n_max = 1,
               skip = 1,
               tokenizer = tokenizer_tsv()) - 6
col_names <- c("chrom", "start", "end","gene", "value", "sign", 1:num_bins)

meta <- read_lines(sense_files[1],n_max = 1)
sample_labels <- str_c(my_names,collapse = "\",\"")
# merge data
tablefile <- NULL
for(i in seq_along(sense_files)){
  if(is.null(tablefile)){
    tablefile <- suppressMessages(read_tsv(
      sense_files[i],
      comment = "@",
      col_names = col_names)) 
  } else {
    tablefile <- suppressMessages(read_tsv(
      sense_files[i],
      comment = "@",
      col_names = col_names)) %>% 
      full_join(tablefile, ., by=c("chrom", "start", "end","gene", "value", "sign"))
  }
  
}
if(length(asense_files)>0){
  asample_labels <- str_c(my_names,"_antisense",collapse = "\",\"")
  sample_labels <- str_c(sample_labels,"\",\"",asample_labels,collapse = "")
  for(i in seq_along(asense_files)){
    tablefile <- suppressMessages(read_tsv(
        asense_files[i],
        comment = "@",
        col_names = col_names)) %>% 
        full_join(tablefile, ., by=c("chrom", "start", "end","gene", "value", "sign"))
    }
}

tablefile[is.na(tablefile)] <- 0
tablefile <- tablefile 
n <- length(sense_files) + length(asense_files)

# fix header

# update sample labels
replacement_string <- paste0("\"sample_labels\":[\"", sample_labels, "\"]")
replacement_group <- paste0("\"group_labels\":[\"", groulab, "\"]")
# Apply gsub with the modified replacement string
meta <- gsub("\"sample_labels\":\\[\"(.*?)\"\\]", replacement_string, meta)
meta <- gsub("\"group_labels\":\\[\"(.*?)\"\\]", replacement_group, meta)
# repeate all numbers n times
meta <- gsub("\\[(\\d+)\\]", paste0("\\[",paste0(rep("\\1",n),collapse = ","),"\\]"), meta)
# repeat ref point n times
meta <- gsub("ref point\":\\[(.*?)\\]", paste0("ref point\":\\[",paste0(rep("\\1",n),collapse = ","),"\\]"), meta)
# update number of genes
meta <- gsub("group_boundaries\":\\[0,\\d+\\]", paste0("group_boundaries\":\\[0,",nrow(tablefile),"\\]"), meta)
# set sample boundaries
num_bins <- str_c(seq(from=num_bins,to = num_bins*n,by = num_bins),collapse = ",")
meta <- gsub("sample_boundaries\":\\[0,\\d+\\]", paste0("sample_boundaries\":\\[0,",num_bins,"\\]"), meta)

write_lines(meta,snakemake@output[[1]])
write_tsv(tablefile, snakemake@output[[1]], col_names = F,append = T)

