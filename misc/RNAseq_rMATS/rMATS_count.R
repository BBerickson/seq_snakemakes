# benjamin.erickson@cuanschutz.edu

# Counting total numbers from rMATS output

#PSI values
# PSI (Percent Spliced In), range from 0 (which would indicate that the exon is never included) to 1 (which would indicate that the exon is always included). 

# * **ID** A unique identifier for this event in this table. Useful for only that reason.
# * **chr** chromosome
# * **strand** strand that the transcript is on (+ or -)
# * **exonStart_0base** the coordinate of the beginning of the alternative exon (using 0-based coordinates)
# * **exonEnd** the coordinate of the end of the alternative exon
# * **upstreamES** the coordinate of the beginning of the exon immediately upstream of the alternative exon
# * **upstreamEE** the coordinate of the end of the exon immediately upstream of the alternative exon
# * **downstreamES** the coordinate of the beginning of the exon immediately downstream of the alternative exon
# * **downstreamEE** the coordinate of the end of the exon immediately downstream of the alternative exon
# * **IJC_SAMPLE_X** the number of read counts that support inclusion of the exon is sample X (four numbers, one for each replicate, each separated by a comma)
# * **SJC_SAMPLE_X** same thing, but for read counts that support the exclusion of the exon
# * **PValue** The pvalue asking if the PSI values for this event between the two conditions is statistically significantly different
# * **FDR** The p value, after it has been corrected for multiple hypothesis testing. This is the significance value you would want to filter on.
# * **IncLevel1** PSI values for the replicates in condition 1 (in this case, condition 1 is RBFOX shRNA).
# * **IncLevel2** PSI values for the replicates in condition 2 (in this case, condition 2 is Control shRNA).
# * **IncLevelDifference** Difference in PSI values between conditions (Condition 1 - Condition 2).

library(tidyverse)
setwd("PROJ/rmats") # set working dir

sigcounts <- function(comp = "rMATS_folder", type = "SE", min_count = 2, max_FDR=0.05, max_diffrence = 0.2, savefiles=FALSE){
  psis <- read_tsv(paste0(comp,'/',type,'.MATS.JC.txt'),show_col_types = FALSE) 
  
  out_length <- length(str_split(psis[1,"IJC_SAMPLE_1"],",")[[1]])
  
  sumofcol <- function(df, col1,col2, colnum) {
    mutate(df, !!{{colnum}} := !! {{col1}} + !! {{col2}})
  }
  
  psis <- psis %>%
    #Split the replicate read counts that are separated by commas into different columns
    separate(., col = IJC_SAMPLE_1, into = paste0('IJC_S1R', 1:out_length), sep = ',', remove = T, convert = T) %>%
    separate(., col = SJC_SAMPLE_1, into = paste0('SJC_S1R', 1:out_length), sep = ',', remove = T, convert = T) %>%
    separate(., col = IJC_SAMPLE_2, into = paste0('IJC_S2R', 1:out_length), sep = ',', remove = T, convert = T) %>%
    separate(., col = SJC_SAMPLE_2, into = paste0('SJC_S2R', 1:out_length), sep = ',', remove = T, convert = T)
  
  # Now sum to get reads in each condition for each event and filter.
  # adding the counts in each sample of inclusion and exclusion 
  psis.filtered <- psis
  for(i in 1:out_length){
    psis.filtered <- sumofcol(psis.filtered,
                              as.name(paste0("IJC_S1R",i)),as.name(paste0("SJC_S1R",i)),
                              as.name(paste0("S1R",i,"counts"))) %>%
      dplyr::filter(.,!!as.name(paste0("S1R",i,"counts")) >= min_count) %>% 
      sumofcol(.,
               as.name(paste0("IJC_S2R",i)),as.name(paste0("SJC_S2R",i)),
               as.name(paste0("S2R",i,"counts"))) %>% 
      dplyr::filter(.,!!as.name(paste0("S2R",i,"counts")) >= min_count)
    
  }
  
  # Defining sensitive exons #only those whose PSI decreases < max_diffrence
  psis.sensitive1 <- dplyr::filter(psis.filtered, FDR < max_FDR,  IncLevelDifference < -max_diffrence) 
  # Defining sensitive exons #only those whose PSI increase > max_diffrence
  psis.sensitive2 <- dplyr::filter(psis.filtered, FDR < max_FDR,  IncLevelDifference > max_diffrence) 
  
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
    less <- "more_inclusion"
    more <- "more_skipping"
  } else if(type == "RI"){
    less <- "more_retained_introns"
    more <- "less_retained_introns"
  } else if(type == "A5SS"){
    less <- "more_A5SS"
    more <- "less_A5SS"
  } else if(type == "A3SS"){
    less <- "more_A3SS"
    more <- "less_A3SS"
  } else if(type == "MXE"){
    less <- "more_MXE"
    more <- "less_MXE"
  }
  if(savefiles){
    psis.sensitive1 %>%  write_tsv(.,paste0(comp,'_',less,'_min_count',min_count,'.MATS.JC.filter.txt'),col_names = T)
    psis.sensitive2 %>%  write_tsv(.,paste0(comp,'_',more,'_min_count',min_count,'.MATS.JC.filter.txt'),col_names = T)
  }
 
  v1 = paste("less",type,comp,sep="_")
  v2 = paste("more",type,comp,sep="_")
  tibble(!!v1 := nrow(psis.sensitive1), !!v2 :=nrow(psis.sensitive2))
  
}

comps <- c("SE","RI", "MXE","A5SS","A3SS")

SE <- NULL
for(i in c(comps)){
  SE <- bind_cols(SE,sigcounts("treatment:control",type = i,min_count = 2, max_FDR=0.05, max_diffrence = 0.2)) # set treatment:control folder name
}

SE
library(ggpubr)
df2 <- pivot_longer(SE,cols=names(SE), names_to = "sample",values_to = "count")
df2 <- df2 %>% mutate(type=str_sub(sample, start = 6),category=str_sub(sample,end = 4)) %>% 
  mutate(type=str_replace(type,"_","&")) %>% 
  separate(type,into = c("type","temp"),sep = "&") %>% 
  select(-temp)

ggbarplot(data = df2, 
          x = "type", y = "count", 
          fill = "category",  # Use 'category' to stack
          color = "category", 
          palette = "jco",  # Predefined color palette
          label = TRUE, 
          position = position_stack()) +  # Stacked bar
  rotate_x_text(angle = 45) +  # Rotate y-axis labels if needed
  labs(x = "Event Type", y = "Count", title = "treatment:control, min_count = 2, max_FDR=0.05, max_diffrence = 0.2")


