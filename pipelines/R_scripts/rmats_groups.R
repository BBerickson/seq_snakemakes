#!/usr/bin/env Rscript

library("tidyverse")

mydir <- snakemake@params[["project"]]
GRPS_UNIQ <- snakemake@params[["mygroup"]] %>% str_split_fixed(.,":",2) %>% c() %>% paste0("-",.,"$")
SAM_GRPS <- snakemake@params[["mysamps"]]
COLS_DICT <- snakemake@params[["mycols"]]

samp1 <- SAM_GRPS[str_detect(SAM_GRPS,GRPS_UNIQ[1])] %>% str_remove(.,GRPS_UNIQ[1]) 
grps <- samp1 %>% paste0(.,"*.bam$")
bams <- list.files(path = mydir, pattern = glob2rx(paste0("^", grps,collapse = "|"))) %>% paste0(mydir,"/",.)
gglist <- tibble(sample=samp1,files=bams,group=str_remove(GRPS_UNIQ[1],"-") %>% str_remove(.,"\\$"))

# Combine pairs with commas and save file
write_file(paste(bams, collapse = ","), snakemake@output[[1]])

samp2 <- SAM_GRPS[str_detect(SAM_GRPS,GRPS_UNIQ[2])] %>% str_remove(.,GRPS_UNIQ[2]) 
grps <- samp2 %>% paste0(.,"*.bam$")
bams <- list.files(path = mydir, pattern = glob2rx(paste0("^", grps,collapse = "|"))) %>% paste0(mydir,"/",.)
gglist <- bind_rows(tibble(sample=samp2,files=bams,group=str_remove(GRPS_UNIQ[2],"-") %>% str_remove(.,"\\$")),
                    gglist)
mycols <- tibble(
  sample = names(COLS_DICT),
  color = map_chr(COLS_DICT,1)
) %>% inner_join(gglist,.,by="sample") %>% distinct(color)

# Convert to hex
mycols <- mycols %>%
  mutate(color = map_chr(color, function(x) {
    if (str_detect(x, "^#")) {
      # Already a hex code
      x
    } else if (str_detect(x, "^[0-9]+,[0-9]+,[0-9]+$")) {
      # RGB string
      rgb_vals <- as.numeric(str_split(x, ",")[[1]])
      rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], maxColorValue = 255)
    } else {
      # Named color
      col2rgb(x) %>%
        as.numeric() %>%
        { rgb(.[1], .[2], .[3], maxColorValue = 255) }
    }
  }))

# Combine pairs with commas and save file
write_file(paste(bams, collapse = ","), snakemake@output[[2]])
write_tsv(gglist,snakemake@output[[3]],col_names = F)
write_tsv(mycols,snakemake@output[[4]],col_names = F)
