# 2025-02-05
setwd("~/Desktop/URLs/220708_U170K")
# load cleavage ULR's in to bentools v9 of the 3' upstream
# 3' table + 1kb:15kb 100bp bins = 25bins
db <- LIST_DATA$table_file %>% select(gene,score,bin,set) %>% 
  mutate(set=str_remove(set,"rep1_|rep2_|rep3_")) %>% 
  mutate(set=str_remove(set,"norm_cleavage_BL_")) %>% 
  mutate(set=str_remove(set,"_sense")) 
# upstream filtering for >= 10 -1kb:0.5kb
upstream <- db %>% 
  filter(bin <= 5) %>% group_by(gene,set) %>% 
  summarise(score=mean(score,na.rm = T),.groups = "drop") %>% 
  filter(score >=10) 
# downstream, filtering for >= 5 0.5kb:1kb
downstream <- db %>% 
  filter(bin > 15,bin <= 20) %>% group_by(gene,set) %>% 
  summarise(score=mean(score,na.rm = T),.groups = "drop") %>% 
  filter(score >=5) 
# inner join temp tables, calculate clevage scores down/up, and filter > 0.2, and make gene list for nascent gene filter
cleavage <- inner_join(upstream,downstream,by=c("gene","set"),suffix = c(".up",".down"))
cleavage %>% mutate(cleavage=score.down/score.up) %>% filter(cleavage > 0.2) %>% group_by(gene) %>%
  mutate(min_value = min(cleavage)) %>%
  ungroup() %>%
  arrange(desc(min_value), gene, cleavage) %>%
  select(-min_value) %>% 
  write_tsv(.,"220708_U170K_Bruseq_cleavage_filter.tsv",col_names = F)
cleavage %>% mutate(cleavage=score.down/score.up) %>% filter(cleavage > 0.2) %>% group_by(gene) %>%
  mutate(min_value = min(cleavage)) %>%
  ungroup() %>%
  arrange(desc(min_value), gene, cleavage) %>%
  select(-min_value) %>% 
  write_tsv(.,"220708_U170K_Bruseq_cleavage_filter_lessthan0.2.tsv",col_names = F)

cleavage %>% mutate(cleavage=score.down/score.up) %>% group_by(gene) %>%
  mutate(min_value = min(cleavage)) %>%
  ungroup() %>%
  arrange(desc(min_value), gene, cleavage) %>%
  select(-min_value) %>% 
  write_tsv(.,"220708_U170K_Bruseq_cleavage_all.tsv")
