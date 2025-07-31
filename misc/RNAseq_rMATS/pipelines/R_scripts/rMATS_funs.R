# benjamin.erickson@cuanschutz.edu



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

# Counting total numbers from rMATS output
sigcounts <- function(PROJ, tc, type = "SE", min_count = 2, max_FDR = 0.05, max_IncDifference = 0.2, savefiles=FALSE){
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

# rMATS outputs
sigoutputs <- function(PROJ,tc, PSI_up = TRUE, type = "SE", min_count = 2, max_FDR = 0.05, max_IncDifference = 0.2, savefiles=FALSE){
  psis <- read_tsv(paste0(PROJ,"/rmats/",tc,'/',type,'.MATS.JC.txt'),show_col_types = FALSE) 
  
  out_length_1 <- length(str_split(psis[1,"IJC_SAMPLE_1"],",")[[1]])
  out_length_2 <- length(str_split(psis[1,"IJC_SAMPLE_2"],",")[[1]])
  
  sumofcol <- function(df, col1,col2, colnum) {
    mutate(df, !!{{colnum}} := !! {{col1}} + !! {{col2}})
  }
  
  psis <- psis %>%
    #Split the replicate read counts that are separated by commas into different columns
    separate(., col = IJC_SAMPLE_1, into = paste0('IJC_S1R', 1:out_length_1), sep = ',', remove = T, convert = T) %>%
    separate(., col = SJC_SAMPLE_1, into = paste0('SJC_S1R', 1:out_length_1), sep = ',', remove = T, convert = T) %>%
    separate(., col = IJC_SAMPLE_2, into = paste0('IJC_S2R', 1:out_length_2), sep = ',', remove = T, convert = T) %>%
    separate(., col = SJC_SAMPLE_2, into = paste0('SJC_S2R', 1:out_length_2), sep = ',', remove = T, convert = T)
  
  # Now sum to get reads in each condition for each event and filter.
  # adding the counts in each sample of inclusion and exclusion 
  psis.filtered <- psis
  for(i in 1:out_length_1){
    psis.filtered <- sumofcol(psis.filtered,
                              as.name(paste0("IJC_S1R",i)),as.name(paste0("SJC_S1R",i)),
                              as.name(paste0("S1R",i,"counts"))) %>%
      dplyr::filter(.,!!as.name(paste0("S1R",i,"counts")) >= min_count)
  }
  for(i in 1:out_length_2){
    psis.filtered <- sumofcol(psis.filtered,
                              as.name(paste0("IJC_S2R",i)),as.name(paste0("SJC_S2R",i)),
                              as.name(paste0("S2R",i,"counts"))) %>% 
      dplyr::filter(.,!!as.name(paste0("S2R",i,"counts")) >= min_count)
    
  }
  
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
  if(savefiles){
    psis.sensitive1 %>%  write_tsv(.,paste0(PROJ,"/rmats/",tc,'_',less,'_min_count',min_count,'.MATS.JC.filter.txt'),col_names = T)
    psis.sensitive2 %>%  write_tsv(.,paste0(PROJ,"/rmats/",tc,'_',more,'_min_count',min_count,'.MATS.JC.filter.txt'),col_names = T)
  }
  
  if(PSI_up){
    return(
      psis.sensitive2
    )
  }else{
    return(
      psis.sensitive1
    )
  }
  
  
  
}

# make gene lists for ggsashimi

makeSashimi <- function(PROJ, tc, type="SE", num_output = 20, min_count = 2, max_FDR=0.05, max_IncDifference = 0.2, gene_Symbol = NULL){
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
  psis <- read_tsv(paste0(PROJ,"/rmats/",tc,'/',type,'.MATS.JC.txt'),show_col_types = FALSE) 
  
  if(type == "A5SS" | type == "A3SS"){
    psis <- psis %>% dplyr::rename(upstreamES=longExonStart_0base,downstreamEE=longExonEnd)
  }
  
  if(!is.null(gene_Symbol)){
    for(i in gene_Symbol){
      psis.filtered <- psis %>% 
        dplyr::filter(geneSymbol == i) %>% 
        group_by(`ID...1`) %>% 
        mutate(my_count = rowSums(across(ends_with("counts")))) %>%
        arrange(IncLevelDifference,desc(my_count),FDR) %>% 
        group_by(geneSymbol) %>% 
        mutate(upstreamES=min(upstreamES),downstreamEE=max(downstreamEE)) %>% 
        ungroup() %>% 
        distinct(upstreamES,downstreamEE,.keep_all = T) %>% 
        slice_head(., n=num_output)
      
      psis.sensitive <- psis.filtered %>%
        #Get rid of columns we aren't really going to use.
        select("geneSymbol","chr","upstreamES","downstreamEE","FDR", "strand") %>%
        mutate(name=paste0(chr,":",upstreamES,"-",downstreamEE)) %>% 
        mutate(geneSymbol=if_else(is.na(geneSymbol),name,paste0(geneSymbol,"-",name)))
      
      psis.sensitive %>% select(name,geneSymbol) %>% write_tsv(.,paste0(PROJ,"/rmats/",tc,"_",i,".txt"),col_names = F)
    }
  } else {
    out_length_1 <- length(str_split(psis[1,"IJC_SAMPLE_1"],",")[[1]])
    out_length_2 <- length(str_split(psis[1,"IJC_SAMPLE_2"],",")[[1]])
    
    sumofcol <- function(df, col1,col2, colnum) {
      mutate(df, !!{{colnum}} := !! {{col1}} + !! {{col2}})
    }
    
    psis <- psis %>%
      #Split the replicate read counts that are separated by commas into different columns
      separate(., col = IJC_SAMPLE_1, into = paste0('IJC_S1R', 1:out_length_1), sep = ',', remove = T, convert = T) %>%
      separate(., col = SJC_SAMPLE_1, into = paste0('SJC_S1R', 1:out_length_1), sep = ',', remove = T, convert = T) %>%
      separate(., col = IJC_SAMPLE_2, into = paste0('IJC_S2R', 1:out_length_2), sep = ',', remove = T, convert = T) %>%
      separate(., col = SJC_SAMPLE_2, into = paste0('SJC_S2R', 1:out_length_2), sep = ',', remove = T, convert = T)
    
    # Now sum to get reads in each condition for each event and filter.
    # adding the counts in each sample of inclusion and exclusion 
    psis.filtered <- psis
    for(i in 1:out_length_1){
      psis.filtered <- sumofcol(psis.filtered,
                                as.name(paste0("IJC_S1R",i)),as.name(paste0("SJC_S1R",i)),
                                as.name(paste0("S1R",i,"counts"))) %>%
        dplyr::filter(.,!!as.name(paste0("S1R",i,"counts")) >= min_count)
    }
    for(i in 1:out_length_2){
      psis.filtered <- sumofcol(psis.filtered,
                                as.name(paste0("IJC_S2R",i)),as.name(paste0("SJC_S2R",i)),
                                as.name(paste0("S2R",i,"counts"))) %>% 
        dplyr::filter(.,!!as.name(paste0("S2R",i,"counts")) >= min_count)
      
    }
    
    # Defining sensitive exons #only those whose PSI decreases < 0
    psis.sensitive <- filter(psis.filtered, FDR < max_FDR, IncLevelDifference < -max_IncDifference) %>% 
      group_by(`ID...1`) %>% 
      mutate(my_count = rowSums(across(ends_with("counts")))) %>%
      arrange(IncLevelDifference,desc(my_count),FDR) %>% 
      group_by(geneSymbol) %>% 
      mutate(upstreamES=min(upstreamES),downstreamEE=max(downstreamEE)) %>% 
      ungroup() %>% 
      distinct(upstreamES,downstreamEE,.keep_all = T) %>% 
      slice_head(., n=num_output)
    
    psis.sensitive <- psis.sensitive %>%
      #Get rid of columns we aren't really going to use.
      select("geneSymbol","chr","upstreamES","downstreamEE","FDR", "strand") %>%
      mutate(name=paste0(chr,":",upstreamES,"-",downstreamEE)) %>% 
      mutate(geneSymbol=if_else(is.na(geneSymbol),name,paste0(geneSymbol,"-",name)))
    
    psis.sensitive %>% filter(strand == "+") %>% select(name,geneSymbol) %>% write_tsv(.,paste0(PROJ,"/rmats/",tc,"_",less,"_pos.txt"),col_names = F)
    psis.sensitive %>% filter(strand == "-") %>% select(name,geneSymbol) %>% write_tsv(.,paste0(PROJ,"/rmats/",tc,"_",less,"_neg.txt"),col_names = F)
    
    psis.sensitive <- filter(psis.filtered, FDR < max_FDR, IncLevelDifference > max_IncDifference) %>% 
      mutate(my_count = rowSums(across(ends_with("counts")))) %>% 
      arrange(desc(IncLevelDifference),desc(my_count),FDR) %>% 
      group_by(geneSymbol) %>% 
      mutate(upstreamES=min(upstreamES),downstreamEE=max(downstreamEE)) %>% 
      ungroup() %>% 
      distinct(upstreamES,downstreamEE,.keep_all = T) %>% 
      slice_head(., n=num_output)
    
    psis.sensitive <- psis.sensitive %>%
      #Get rid of columns we aren't really going to use.
      select("geneSymbol","chr","upstreamES","downstreamEE","FDR", "strand") %>%
      mutate(name=paste0(chr,":",upstreamES,"-",downstreamEE)) %>% 
      mutate(geneSymbol=if_else(is.na(geneSymbol),name,paste0(geneSymbol,"-",name)))
    
    psis.sensitive %>% filter(strand == "+") %>% select(name,geneSymbol) %>% write_tsv(.,paste0(PROJ,"/rmats/",tc,"_",more,"_pos.txt"),col_names = F)
    psis.sensitive %>% filter(strand == "-") %>% select(name,geneSymbol) %>% write_tsv(.,paste0(PROJ,"/rmats/",tc,"_",more,"_neg.txt"),col_names = F)
  }
}



