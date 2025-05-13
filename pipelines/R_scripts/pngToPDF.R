#!/usr/bin/env Rscript

library("png")
library("grid")
library("gridExtra")

# Read input and parameters from Snakemake
pngFiles <- snakemake@input
chunkSize <- snakemake@params[["chunkSize"]]
nsize <- length(pngFiles)

# Determine the number of columns for the grid layout
mycols <- ifelse(chunkSize == 1, 1, 2)

# Open the PDF device
pdf(snakemake@output[[1]])

# Loop through the chunks and create the grid layout
for(i in seq(1, nsize, by = chunkSize)){  
  # Read the PNG files for the current chunk
  rl <- lapply(pngFiles[i:min(nsize, (i + chunkSize - 1))], png::readPNG)
  # Convert the PNG files to raster grobs
  gl <- lapply(rl, grid::rasterGrob)
  # Arrange the grobs in a grid layout
  gridExtra::grid.arrange(grobs = gl, ncol = mycols)
}

# Close the PDF device
dev.off()



