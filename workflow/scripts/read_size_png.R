#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library("tidyverse")

samp_file <- args[1]
samp_png <- args[2]

df <- read_tsv(samp_file,skip = 9) %>% dplyr::rename(Length=`#Length`) %>% filter(Length < 1001)

# Plot the read length histogram
p <- ggplot(df, aes(x=Length, y=reads)) +
  geom_bar(stat="identity", fill="steelblue", color="black") +
  labs(title="Read Length Distribution", x="Read Length", y="Number of Reads") +
  theme_minimal()

png(samp_png, width = 700, height = 525)
print(p)  # Print the ggplot object to the device
dev.off()  # Close the device

