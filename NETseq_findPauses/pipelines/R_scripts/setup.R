

# Set up consents
library(here)
library(rlang)
library(vroom)
library(scales)
library(colorspace)
library(cowplot)
# library(ggseqlogo)       # To create sequence logos
library(valr)            # To create sequence logos
# library(bedr)            # To create sequence logos

# Gene lengths for filtering
MIN_GENE_LEN <- 2
MIN_ZONE_GENE_LEN <- 4

# "Theme"
txt_mplyr <- 1.5

ttl_pt <- 10 * txt_mplyr
txt_pt <- 8  * txt_mplyr
ln_pt  <- 0.5

theme_info <- theme_cowplot() +
  theme(
    plot.title       = element_text(face = "plain"),
    strip.background = element_rect(fill = NA),
    legend.title     = element_text(size = ttl_pt),
    legend.text      = element_text(size = ttl_pt),
    strip.text       = element_text(size = ttl_pt),
    axis.title       = element_text(size = ttl_pt),
    axis.text        = element_text(size = txt_pt),
    axis.line        = element_blank(),
    panel.border     = element_rect(fill = NA, color = "black", linewidth = ln_pt),
    axis.ticks       = element_line(color = "black", linewidth = ln_pt)
  )

theme_colors <- c(
  "#225ea8",  # blue
  "#e31a1c",  # red
  "#238443",  # green
  "#ec7014",  # orange
  "#8c6bb1",  # purple
  "#41b6c4",  # aqua
  "#737373"   # grey
)

# Columns for bed files
bed_cols <- c(
  "chrom", "start", 
  "end",   "name", 
  "score", "strand"
)

box_theme <- theme_info +
  theme(
    plot.title       = element_text(size = ttl_pt),
    legend.position  = "bottom",
    legend.direction = "vertical",
    legend.title     = element_blank(),
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.line        = element_blank(),
    panel.border     = element_rect(fill = NA, color = "black")
  )
