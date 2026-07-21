#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)


PROJ <- args[1] # same as PROJ: in samples.yaml
tc <- args[2] # "treatment:control" same as GROUPS_COMP: in samples_rAMTS.yaml
min_count <- as.numeric(args[3]) # >= min_count average reads per replicate (inclusion + skipping) in at least one condition.
max_FDR <- as.numeric(args[4]) # Standard 5% significance threshold for FDR adjusted pValue
max_IncDifference <- as.numeric(args[5]) # PSI change threshold
comps <- unlist(str_split(args[6], "\\s+")) # c("SE","RI", "MXE","A5SS","A3SS") # splicing types to loop over
num_plots <- as.numeric(args[7])
pdf_out <- args[8]

sigcounts <- function(PROJ, tc, type = "SE", min_count = 20, max_FDR = 0.05, max_IncDifference = 0.05, savefiles=TRUE){
  psis <- read_tsv(paste0(PROJ,"/rmats/",tc,'/',type,'.MATS.JC.txt'),show_col_types = FALSE) 
  
  out_length_1 <- length(str_split(psis[1,"IJC_SAMPLE_1"],",")[[1]])
  out_length_2 <- length(str_split(psis[1,"IJC_SAMPLE_2"],",")[[1]])
  
  # Calculate scaled thresholds based on number of replicates
  min_count_condition1 <- min_count * out_length_1
  min_count_condition2 <- min_count * out_length_2
  
  psis <- psis %>%
    #Split the replicate read counts that are separated by commas into different columns
    separate(., col = IJC_SAMPLE_1, into = paste0('IJC_S1R', 1:out_length_1), sep = ',', remove = T, convert = T) %>%
    separate(., col = SJC_SAMPLE_1, into = paste0('SJC_S1R', 1:out_length_1), sep = ',', remove = T, convert = T) %>%
    separate(., col = IJC_SAMPLE_2, into = paste0('IJC_S2R', 1:out_length_2), sep = ',', remove = T, convert = T) %>%
    separate(., col = SJC_SAMPLE_2, into = paste0('SJC_S2R', 1:out_length_2), sep = ',', remove = T, convert = T)
  
  # Calculate total counts per condition and filter
  # Sum IJC and SJC counts across all replicates for each condition
  ijc_s1_cols <- paste0('IJC_S1R', 1:out_length_1)
  sjc_s1_cols <- paste0('SJC_S1R', 1:out_length_1)
  ijc_s2_cols <- paste0('IJC_S2R', 1:out_length_2)
  sjc_s2_cols <- paste0('SJC_S2R', 1:out_length_2)
  
  psis.filtered <- psis %>%
    rowwise() %>%
    mutate(
      # Total counts per condition (inclusion + exclusion)
      S1_total_counts = sum(c_across(all_of(c(ijc_s1_cols, sjc_s1_cols))), na.rm = TRUE),
      S2_total_counts = sum(c_across(all_of(c(ijc_s2_cols, sjc_s2_cols))), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # Keep events where at least one condition passes the scaled threshold
    filter(S1_total_counts >= min_count_condition1 | S2_total_counts >= min_count_condition2)
  
  # Defining sensitive exons #only those whose PSI decreases < max_IncDifference
  psis.sensitive1 <- dplyr::filter(psis.filtered, FDR < max_FDR,  IncLevelDifference < -max_IncDifference) 
  # Defining sensitive exons #only those whose PSI increase > max_IncDifference
  psis.sensitive2 <- dplyr::filter(psis.filtered, FDR < max_FDR,  IncLevelDifference > max_IncDifference) 
  
  # IncLevelDifference = (Condition 1 - Condition 2).
  # IncLevelDifference > 0, explanation for the different PSI types :
  # SE less exon skiping in condition 1
  # MXE:
  # + include exon1 skip exon 2
  # - include exon2 skip exon 1
  # A5SS more downstream
  # A3SS more upstream
  # RI retains intron
  
  if(type == "SE"){
    less <- "more_skipping" # IncLevelDifference < -threshold
    more <- "more_inclusion"  # IncLevelDifference > threshold
  } else if(type == "RI"){
    less <- "less_retained_introns"
    more <- "more_retained_introns"
  } else if(type == "A5SS"){
    less <- "more_downstream_A5SS"
    more <- "more_upstream_A5SS"
  } else if(type == "A3SS"){
    less <- "more_downstream_A3SS"
    more <- "more_upstream_A3SS"
  } else if(type == "MXE"){
    less <- "prefer_exon2"
    more <- "prefer_exon1"
  }
  
  if(type == "SE"){
    se.sensitive <- psis.sensitive1 %>%
      group_by(`ID...1`) %>% 
      mutate(my_count = rowSums(across(ends_with("counts")))) %>%
      arrange(IncLevelDifference,desc(my_count),FDR) %>% 
      group_by(geneSymbol) %>% 
      mutate(upstreamES=min(upstreamES),downstreamEE=max(downstreamEE)) %>% 
      ungroup() %>% 
      distinct(upstreamES,downstreamEE,.keep_all = T) %>% 
      slice_head(., n=num_plots)
    
    se.sensitive <- se.sensitive %>%
      #Get rid of columns we aren't really going to use.
      select("geneSymbol","chr","upstreamES","downstreamEE","FDR", "strand") %>%
      mutate(name=paste0(chr,":",upstreamES,"-",downstreamEE)) %>% 
      mutate(geneSymbol=if_else(is.na(geneSymbol),name,paste0(geneSymbol,"-",name)))
    
    se.sensitive %>% filter(strand == "+") %>% select(name,geneSymbol) %>% write_tsv(.,paste0(PROJ,"/rmats/",tc,"_",less,"_pos.txt"),col_names = F)
    se.sensitive %>% filter(strand == "-") %>% select(name,geneSymbol) %>% write_tsv(.,paste0(PROJ,"/rmats/",tc,"_",less,"_neg.txt"),col_names = F)
    
    se.sensitive <- psis.sensitive2 %>%
      group_by(`ID...1`) %>% 
      mutate(my_count = rowSums(across(ends_with("counts")))) %>%
      arrange(IncLevelDifference,desc(my_count),FDR) %>% 
      group_by(geneSymbol) %>% 
      mutate(upstreamES=min(upstreamES),downstreamEE=max(downstreamEE)) %>% 
      ungroup() %>% 
      distinct(upstreamES,downstreamEE,.keep_all = T) %>% 
      slice_head(., n=num_plots)
    
    se.sensitive <- se.sensitive %>%
      #Get rid of columns we aren't really going to use.
      select("geneSymbol","chr","upstreamES","downstreamEE","FDR", "strand") %>%
      mutate(name=paste0(chr,":",upstreamES,"-",downstreamEE)) %>% 
      mutate(geneSymbol=if_else(is.na(geneSymbol),name,paste0(geneSymbol,"-",name)))
    
    se.sensitive %>% filter(strand == "+") %>% select(name,geneSymbol) %>% write_tsv(.,paste0(PROJ,"/rmats/",tc,"_",more,"_pos.txt"),col_names = F)
    se.sensitive %>% filter(strand == "-") %>% select(name,geneSymbol) %>% write_tsv(.,paste0(PROJ,"/rmats/",tc,"_",more,"_neg.txt"),col_names = F)
  }
  
  if(savefiles){
    psis.sensitive1 %>%  write_tsv(.,paste0(PROJ,"/rmats/",tc,'_',less,'_min_count',min_count,'.MATS.JC.filter.txt'),col_names = T)
    psis.sensitive2 %>%  write_tsv(.,paste0(PROJ,"/rmats/",tc,'_',more,'_min_count',min_count,'.MATS.JC.filter.txt'),col_names = T)
  }
  
  # Return a tibble with the desired structure
  total_events <- nrow(psis)
  filtered_events <- sum(nrow(psis.sensitive2), nrow(psis.sensitive1))
  percent_passed <- round((filtered_events / total_events) * 100, 2)
  
  return(
    tibble(
      sample = tc,
      count = c(nrow(psis.sensitive2), nrow(psis.sensitive1)),
      type = type,
      category = c(more, less),
      direction = c("more", "less"),
      total_events = total_events,
      filtered_events = filtered_events,
      percent_passed = percent_passed
    )
  )
}

##### bar plots #####
SE <- list()
for(i in c(comps)){
  SE[[i]] <- sigcounts(PROJ,tc,type = i,min_count, max_FDR, max_IncDifference) # set treatment:control folder name
}

# Combined and Define which categories are "less" and should be negative
db <- bind_rows(SE)%>%
  mutate(
    signed_count = ifelse(direction == "less", -count, count)
  )
# set levels for plot order
db$type <- factor(db$type, levels = comps)
db$category <- factor(db$category, levels = db$category)


# Create bar plot
p <- ggplot(db, aes(x = type, y = signed_count, fill = category)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_brewer(palette = "Paired") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_text(aes(label = abs(signed_count)), vjust = ifelse(db$signed_count >= 0, -0.2, 1.2), size = 5) +
  labs(x = "PSI Type", y = bquote(Delta~"PSI count:"~"[" * .(gsub(":", " - ", tc)) * "]"), 
       fill = str_split_fixed(tc,":",2)[1],
       title = str_replace(tc,":"," vs "),
       subtitle = paste0("min ave reads/rep = ",min_count, ", max_FDR = ",max_FDR, ", IncDifference = ",max_IncDifference)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 1),
        axis.line.y = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1))

# Group by type and calculate the average percent_passed
percent_plot_data <- db %>%
  group_by(type) %>%
  summarise(percent_passed = mean(percent_passed))

# Create the plot
ppd <- ggplot(percent_plot_data, aes(x = type, y = percent_passed, fill = type)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +
  geom_text(aes(label = paste0(round(percent_passed, 2), "%")), 
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Percentage of Events Passed Filtering by PSI Type",
    x = "PSI Type",
    y = "Percent Passed (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "none"
  )


table_plot <- tableGrob(db %>% select(-direction,-signed_count,-total_events, -filtered_events) %>% separate(sample,into = c("treatment","condition")) %>% mutate(category=paste(category,"in",treatment)))

pdf(pdf_out, width = 10, height = 7.5)
print(p)
print(ppd)
grid.arrange(table_plot)
dev.off()

