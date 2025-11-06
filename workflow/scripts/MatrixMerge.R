#!/usr/bin/env Rscript
# merge deeptools computeMatrix files 
library("tidyverse")
library("valr")
args = commandArgs(trailingOnly=TRUE)
infile <- strsplit(args[1], " ")[[1]]
outfile <- args[2]
gl <- NULL
nickname <- NULL
groulab <- NULL
num_bins <- count_fields(infile[1], n_max = 1, skip = 1, tokenizer = tokenizer_tsv()) -6
for(i in seq_along(infile)){
  infile_bins <- count_fields(infile[i], n_max = 1, skip = 1, tokenizer = tokenizer_tsv()) -6
  if(infile_bins == num_bins){
    meta <- read_lines(infile[i],n_max = 1)
    
    gl <- suppressMessages(read_tsv(
      infile[i],
      comment = "@",
      col_names = c("chrom", "start", "end","gene", "value", "sign", 1:num_bins))) %>% 
      {if (!is.null(gl)) full_join(gl,.,by=c("chrom", "start", "end","gene", "value", "sign")) else .}
    nickname <- c(nickname, meta %>% str_extract('sample_labels":\\["([^"]+)"', group = 1) %>% 
                    str_split_fixed("_aligned_", 2) %>% .[, 1]) 
    groulab <- c(groulab, meta %>% str_extract('group_labels":\\["([^"]+)"', group = 1))
  }
}
gl[is.na(gl)] <- 0
# gather info for filtering 
mylist <- c("ref point:","upstream:","downstream:","unscaled 5 prime:","unscaled 3 prime:")
mm <- meta %>% str_remove_all("[@{}]|\\]|\\[") %>% str_split(",",simplify = T) %>% 
  str_replace_all(.,fixed('\"'),"")
mm <- mm[str_detect(mm,paste(mylist,collapse = "|"))] 
upstream <- mm[str_detect(mm,"upstream")] %>% str_remove("upstream:") %>% as.numeric()
downstream <- mm[str_detect(mm,"downstream")] %>% str_remove("downstream:") %>% as.numeric()
un5 <- mm[str_detect(mm,"unscaled 5 prime:")] %>% str_remove("unscaled 5 prime:") %>% as.numeric()
un3 <- mm[str_detect(mm,"unscaled 3 prime:")] %>% str_remove("unscaled 3 prime:") %>% as.numeric()
if(any(str_detect(mm,"ref point:TSS"))){
  sep_bins <- max(upstream,1)
  gene_length <- max(downstream,10)
} else if(any(str_detect(mm,"ref point:TES"))) {
  sep_bins <- max(downstream,1)
  gene_length <- max(upstream,10)
} else{
  sep_bins <- max(c(upstream,downstream),1)
  gene_length <- max(sum(c(un5, un3)),10)
}

gl <- gl %>% 
  bed_cluster(max_dist = as.numeric(sep_bins))%>% 
  filter(!(duplicated(.id) | duplicated(.id, fromLast=TRUE))) %>% 
  filter(end-start >= as.numeric(gene_length)) %>% 
  mutate(score=rowMeans(across(contains(".x")))) %>% arrange(score) %>% 
  slice_tail(prop = 0.90) %>% select(-.id,-score) %>% bed_sort() 
  

n <- length(infile) 
# Calculate number of genes
ngenes <- nrow(gl)
# fix header
sample_labels <- str_c(nickname,collapse = "\",\"")
replacement_string <- paste0("\"sample_labels\":[\"", sample_labels, "\"]")
# Add gene count to group label
replacement_group <- paste0("\"group_labels\":[\"", paste("n =", ngenes, unique(groulab)[1]), "\"]")
# update sample labels
# Apply gsub with the modified replacement string
meta <- gsub("\"sample_labels\":\\[\"(.*?)\"\\]", replacement_string, meta)
meta <- gsub("\"group_labels\":\\[\"(.*?)\"\\]", replacement_group, meta)
# repeate all numbers n times
meta <- gsub("\\[(\\d+)\\]", paste0("\\[",paste0(rep("\\1",n),collapse = ","),"\\]"), meta)
# repeat ref point n times
meta <- gsub("ref point\":\\[(.*?)\\]", paste0("ref point\":\\[",paste0(rep("\\1",n),collapse = ","),"\\]"), meta)
# update number of genes
meta <- gsub("group_boundaries\":\\[0,\\d+\\]", paste0("group_boundaries\":\\[0,",nrow(gl),"\\]"), meta)
# set sample boundaries
num_bins <- str_c(seq(from=num_bins,to = num_bins*n,by = num_bins),collapse = ",")
meta <- gsub("sample_boundaries\":\\[0,\\d+\\]", paste0("sample_boundaries\":\\[0,",num_bins,"\\]"), meta)
write_lines(meta,outfile)
write_tsv(gl, outfile, col_names = F,append = T)
