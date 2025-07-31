# benjamin.erickson@cuanschutz.edu

# rMATS results

PROJ="180629_SY351" # same as PROJ: in samples.yaml
tc <- "treatment:control" # same as GROUPS_COMP: in samples_rAMTS.yaml
setwd() 

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(gridExtra)
source("pipelines/R_scripts/rMATS_funs.R")

min_count <- 2 # min_count × number of replicates total reads needed to pass filter
max_FDR <- 0.05 # Standard 5% significance threshold for FDR adjusted pValue
max_IncDifference <- 0.2 # PSI change threshold
comps <- c("SE","RI", "MXE","A5SS","A3SS") # splicing types to loop over

##### bar plots #####
SE <- list()
for(i in c(comps)){
  SE[[i]] <- sigcounts(PROJ,tc,type = i,min_count, max_FDR, max_IncDifference, savefiles=FALSE) # set treatment:control folder name
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
ggplot(db, aes(x = type, y = signed_count, fill = category)) +
  geom_bar(stat = "identity", position = "identity") +
  scale_fill_brewer(palette = "Paired") +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  geom_text(aes(label = abs(count)), vjust = ifelse(db$count > 0, -0.2, 1.5), size = 5) +
  labs(x = "PSI Type", y = bquote(Delta~"PSI count:"~"[" * .(gsub(":", " - ", tc)) * "]"), 
       fill = str_split_fixed(tc,":",2)[1],
       title = str_replace(tc,":"," vs "),
       subtitle = paste0("min_count = ",min_count, ", max_FDR = ",max_FDR, ", IncDifference = ",max_IncDifference)) +
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
ggplot(percent_plot_data, aes(x = type, y = percent_passed, fill = type)) +
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

# Display table
grid.arrange(table_plot)

#### sashimi ####

# make gene lists for ggsashimi, num_output will give that many for up and down deltaPSI
makeSashimi(PROJ,tc, min_count, max_FDR, max_IncDifference, num_output = 12) 
## optional usage, will make 1 file for each item 
# makeSashimi(PROJ,tc, gene_Symbol = c("MADD","KLHL9")) 
