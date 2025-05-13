#!/usr/bin/env Rscript

library("tidyverse")
library("stringi")
library("DESeq2")
library("ggplot2")
#my_files <- list.files(pattern = "hg38_featureCounts.tsv$")

my_files <- snakemake@input[["FC"]]
# file names and groupings
mydict <- snakemake@params[["groupkey"]]
# mydict <-  list('PNUTS-ctd'=c('PP1H66KmD_CTD1', 'PP1H66KpD_CTD1'), 'PNUTS-ctd_2'=c('PP1H66KmD_CTD2', 'PP1H66KpD_CTD1'), 'PNUTS-Ser5_1'=c('PP1H66KmD_Ser5_1', 'PP1H66KpD_Ser5_1'), 'PNUTS-Ser5_2'=c('PP1H66KmD_Ser5_2', 'PP1H66KpD_Ser5_2'))

# read in and join data
sub_reads <- read_tsv(my_files[[1]],comment = "#",show_col_types = FALSE)
# set col names

for(i in seq_along(my_files)[-1]) {
  
  sub_reads <- read_tsv(my_files[[i]],comment = "#",show_col_types = FALSE) %>% 
    full_join(sub_reads,., by = names(.)[-length(names(.))])
}

if(any(names(sub_reads) %in% "gene_biotype")){
  sub_reads <- sub_reads %>% dplyr::mutate(Geneid = gene_biotype) %>% dplyr::select(-gene_biotype,-gene_name)
}

# trim down name 
sub_reads <- sub_reads %>% rename_with(~ str_replace(., ".*/", "")) 
# Function to rename columns based on matching patterns
rename_columns <- function(df, dict) {
  dict <- unlist(dict)
  colnames(df) <- sapply(colnames(df), function(col) {
    for (pattern in dict) {
      if (startsWith(col, pattern)) {
        return(pattern)
      }
    }
    return(col)
  })
  return(df)
}
sub_reads <- rename_columns(sub_reads, mydict)

# prepare for DESeq matrix conversion 
sub_reads <- sub_reads %>% 
  dplyr::mutate(gene=paste(Chr,Start,End,Strand,Geneid,sep = ".")) %>% 
  distinct(gene,.keep_all = T) %>% 
  dplyr::select(-Chr,-Start,-End,-Strand,-Geneid,-Length ) 
countData <- sub_reads %>% 
  column_to_rownames(var="gene")

design_formula <- ~ 1

colData <-  stack(mydict) 
colnames(colData) <- c("condition", "group")
row.names(colData) <- colData$condition  
countData <- countData %>% dplyr::select(colData$condition)

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = design_formula) # independent conditions, add group param and update 

# Normalize the data
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)  # Variance stabilizing transformation

# Calculate the PCA and create the plot
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create the PCA plot
p <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  xlim(c(-1, 1) * max(c(abs(pcaData$PC2),abs(pcaData$PC1)))) +  # expand the limits
  ylim(c(-1, 1) * max(c(abs(pcaData$PC2),abs(pcaData$PC1)))) +
  guides(color = guide_legend("Samples"))

#p <- plotPCA(vsd, intgroup=("condition"))

png(snakemake@output[[1]], width = 700, height = 525)
print(p)  # Print the ggplot object to the device
dev.off()  # Close the device

# # save pdf
# pdf(snakemake@output[[1]])
# plotPCA(vsd, intgroup=("condition"))
# dev.off()

