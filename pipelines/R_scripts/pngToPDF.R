#!/usr/bin/env Rscript

library("png")
library("grid")
library("gridExtra")
args = commandArgs(trailingOnly=TRUE)
# Read input and parameters from Snakemake
pngFiles <- strsplit(args[1], " ")[[1]]
chunkSize <- as.numeric(args[2])
title_text <- args[3]  # Add title as 4th argument
outfile <- args[4]


nsize <- length(pngFiles)
# Determine the number of columns for the grid layout
mycols <- ifelse(chunkSize == 1, 1, 2)

# Open the PDF device
pdf(outfile)

# Loop through the chunks and create the grid layout
for(i in seq(1, nsize, by = chunkSize)){  
  # Read the PNG files for the current chunk
  rl <- lapply(pngFiles[i:min(nsize, (i + chunkSize - 1))], png::readPNG)
  # Convert the PNG files to raster grobs
  gl <- lapply(rl, grid::rasterGrob)
  
  # Create title grob
  title_grob <- textGrob(title_text, gp = gpar(fontsize = 12, fontface = "bold"))
  
  # Arrange with title at top
  gridExtra::grid.arrange(title_grob, 
                          arrangeGrob(grobs = gl, ncol = mycols), 
                          ncol = 1, heights = c(1, 10))
}

# Close the PDF device
dev.off()



