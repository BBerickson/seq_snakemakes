#!/usr/bin/env Rscript
library("png")
library("grid")
library("gridExtra")

args <- commandArgs(trailingOnly = TRUE)

# Args: 1=input (dir or space-separated files), 2=chunkSize, 3=title, 4=outfile
input_arg  <- args[1]
chunkSize  <- as.numeric(args[2])
title_text <- args[3]
outfile    <- args[4]

# Resolve input — directory glob or explicit file list
if (dir.exists(input_arg)) {
  message("Input is a directory, globbing PNGs...")
  pngFiles <- sort(list.files(input_arg, pattern = "\\.png$", full.names = TRUE))
} else {
  message("Input is a file list...")
  pngFiles <- strsplit(input_arg, " ")[[1]]
}

# Validate
if (length(pngFiles) == 0) stop("No PNG files found in: ", input_arg)

nsize  <- length(pngFiles)
mycols <- ifelse(chunkSize == 1, 1, 2)

message("Processing ", nsize, " PNGs into: ", outfile)

pdf(outfile)

for (i in seq(1, nsize, by = chunkSize)) {
  chunk <- pngFiles[i:min(nsize, (i + chunkSize - 1))]
  
  rl <- lapply(chunk, png::readPNG)
  gl <- lapply(rl, grid::rasterGrob)
  
  title_grob <- grid::textGrob(
    title_text,
    gp = grid::gpar(fontsize = 12, fontface = "bold")
  )
  
  gridExtra::grid.arrange(
    title_grob,
    gridExtra::arrangeGrob(grobs = gl, ncol = mycols, padding = grid::unit(5, "mm")),
    ncol    = 1,
    heights = c(1, 10)
  )
}

dev.off()
message("Done.")