# 2022-10-20
# make bed file from gene list from bentools
library(tidyverse)
setwd("")
genelist <- read_tsv("genelist.txt",col_names = c("name","value"),comment = "#") %>% filter(name != "gene")

hg38 <- read_bed("/beevol/home/erickson/ref/hg38/ref/hg38_gene.bed",n_fields = 6)
hg38 %>% filter(name %in% genelist$name) %>% 
  write_tsv(.,"genelist.txt.bed",col_names = F)


# 2022-10-12
# make bed file from gene list of top 25% to make filtred matrix files


library("tidyverse")
library("valr")
#bsub < run.sh
setwd("~/Human/mNetSeq/220215_220218_Xrn2_CTD/heatmap")
num_bins <-
  count_fields("snRNA.matrix.gz",
               n_max = 1,
               skip = 1,
               tokenizer = tokenizer_tsv())

tablefile <- suppressMessages(read_tsv(
  "snRNA.matrix.gz",
  comment = "@",
  col_names = c("chrom", "start", "end","gene", "value", "sign", 1:(num_bins - 6)))) 
tablefile2 <- tablefile %>% gather(., bin, score, 7:(all_of(num_bins)))
tablefile2 %>% group_by(chrom,start,end,gene,sign) %>% summarise(score=sum(score),.groups="drop") %>% 
  arrange(desc(score)) %>% slice(1:floor(nrow(.)*.25)) %>% select(chrom,start,end,gene,score,sign) %>%
  filter(sign == "+") %>% 
  write_tsv(.,"snRNA_top25_pos.bed",col_names = F)

tablefile2 %>% group_by(chrom,start,end,gene,sign) %>% summarise(score=sum(score),.groups="drop") %>% 
  arrange(desc(score)) %>% slice(1:floor(nrow(.)*.25)) %>% select(chrom,start,end,gene,score,sign) %>%
  filter(sign == "-") %>% 
  write_tsv(.,"snRNA_top25_neg.bed",col_names = F)

tablefile2 %>% group_by(chrom,start,end,gene,sign) %>% summarise(score=sum(score),.groups="drop") %>% 
  arrange(desc(score)) %>% slice(1:floor(nrow(.)*.25)) %>% select(chrom,start,end,gene,score,sign) %>%
  write_tsv(.,"snRNA_top25.bed",col_names = F)
