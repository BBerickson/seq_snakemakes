#!/usr/bin/env Rscript

library("tidyverse")

df <- read_tsv(snakemake@input[[1]],skip = 9) %>% dplyr::rename(Length=`#Length`) %>% filter(Length < 1001)

# Plot the read length histogram
p <- ggplot(df, aes(x=Length, y=reads)) +
  geom_bar(stat="identity", fill="steelblue", color="black") +
  labs(title="Read Length Distribution", x="Read Length", y="Number of Reads") +
  theme_minimal()

png(snakemake@output[[1]], width = 700, height = 525)
print(p)  # Print the ggplot object to the device
dev.off()  # Close the device

