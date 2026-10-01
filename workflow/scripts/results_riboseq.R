#!/usr/bin/env Rscript
# Forked from results.R -- only the pre-alignment trimming-stats block
# changes (Trim Galore's single-end report format instead of bbduk/clumpify
# or cutadapt), plus an optional UMI-dedup join. Everything else (alignment
# rate, subsample fractions, featureCounts join, fragment-size join) is
# copied as-is.

args <- commandArgs(trailingOnly = TRUE)

library("tidyverse")
library("jsonlite")

# arguments
trim_galore_file <- args[1]  # PROJ + "/stats/" + PROJ + "_trim_galore.tsv"
aligner_file     <- args[2]  # PROJ + "/stats/" + PROJ + "_" + INDEX_MAP + "_aligned.tsv"
dedup_file       <- args[3]  # PROJ + "/stats/" + PROJ + "_" + INDEX_MAP + "_UMI_dedup.tsv", or "" (non-UMI variant)
FC_sum_files     <- strsplit(args[4], " ")[[1]]
FC_files         <- strsplit(args[5], " ")[[1]]

index_map <- args[6]
indexs <- strsplit(args[7], " ")[[1]]
indexs <- indexs[order(nchar(indexs), decreasing = TRUE)] %>% paste0(collapse = "|")
sam_new <- args[8]
mydir <- args[9]
filter_type <- args[10]
outfile <- args[11]

sam_new <- fromJSON(gsub("'", '"', sam_new))
sn <- enframe(sam_new, name = "new_name", value = "sample") %>%
  unnest(sample) %>%
  group_by(new_name) %>%
  mutate(new_name = ifelse(row_number() == 1, new_name, NA)) %>%
  ungroup()

# pre-alignment: Trim Galore's single-end report (wraps cutadapt, but uses SE
# wording -- "Total reads processed"/"Reads written (passing filters)", not
# the PE "Total read pairs processed"/"Pairs written (passing filters)").
pre_alignment <- read_tsv(trim_galore_file, col_names = c("sample", "type", "value"), show_col_types = FALSE) %>%
  dplyr::mutate(type = str_replace_all(type, " ", "_")) %>%
  spread(type, value) %>%
  dplyr::mutate(Total_reads = as.double(Total_reads_processed)) %>%
  dplyr::rename(passing_filters_TrimGalore = `Reads_written_(passing_filters)`) %>%
  dplyr::select(sample, Total_reads, passing_filters_TrimGalore)

if (nzchar(dedup_file)) {
  dedup <- read_tsv(dedup_file, col_names = c("sample", "type", "value"), show_col_types = FALSE) %>%
    dplyr::mutate(sample = str_remove(sample, paste0("_", index_map, "_UMI"))) %>%
    dplyr::mutate(type = str_replace_all(type, " ", "_")) %>%
    dplyr::filter(type == "Input_Reads" | type == "Number_of_reads_out") %>%
    spread(type, value) %>%
    dplyr::mutate(Unique_UMI_reads = paste0(Number_of_reads_out, "(", paste0(round(Number_of_reads_out / Input_Reads * 100, 1), "%"), ")")) %>%
    dplyr::select(sample, Unique_UMI_reads)
  pre_alignment <- full_join(pre_alignment, dedup, by = "sample")
}

# alignment
aligner <- read_tsv(aligner_file, col_names = c("index", "sample", "type", "alignment_rate"), show_col_types = FALSE) %>%
  dplyr::filter(type == "overall_alignment_rate") %>%
  dplyr::select(sample, index, alignment_rate) %>%
  pivot_wider(names_from = index, values_from = alignment_rate, names_glue = "{index}_alignment_rate")

sn <- left_join(sn, pre_alignment, by = "sample") %>%
  left_join(., aligner, by = "sample") %>%
  arrange(new_name) %>%
  mutate(new_name = if_else(is.na(new_name), sample, new_name)) %>%
  select(where(~!all(is.na(.))))

mydirf <- paste0(mydir, "/stats")
mydirc <- paste0(mydir, "/counts")

# feature Counts summary file
if (!is_empty(FC_sum_files)) {
  FC <- NULL
  for (i in seq_along(sn$sample)) {
    FC <- read_tsv(FC_sum_files[str_detect(FC_sum_files, sn$sample[i])],
                   col_names = c("type", "value"), show_col_types = FALSE) %>%
      dplyr::mutate(sample = sn$sample[i]) %>%
      dplyr::filter(!str_detect(type, filter_type)) %>%
      bind_rows(FC)
  }
  FC <- FC %>% distinct(type, sample, .keep_all = TRUE) %>% spread(type, value) %>%
    dplyr::select(sample, distinct(FC, type)$type)
  sn <- left_join(sn, FC, by = "sample") %>% select(where(~!all(is.na(.))))
}

# feature Counts file
if (!is_empty(FC_files)) {
  FCS <- list()
  col_names <- c("ID", "Chr", "Start", "End", "Strand", "Length", "gene_name", "biotype", "count")
  for (i in seq_along(sn$sample)) {
    sub_reads <- read_tsv(FC_files[str_detect(FC_files, sn$sample[i])], comment = "#", show_col_types = FALSE)
    names(sub_reads) <- col_names
    FCS[[i]] <- sub_reads %>% mutate(sample = sn$sample[i]) %>%
      group_by(sample, biotype) %>% summarise(counts = sum(count, na.rm = TRUE), .groups = "drop_last") %>%
      mutate(percentage = (counts / sum(counts, na.rm = TRUE)) * 100) %>%
      mutate(category = case_when((biotype %in% c("protein_coding", "snRNA", "snoRNA", "rRNA")) ~ biotype,
                                   (percentage <= 2) ~ "Other",
                                   (percentage > 2) ~ biotype)) %>%
      group_by(sample, category) %>%
      summarise(percentage = sum(percentage), .groups = "drop") %>%
      arrange(sample)
  }
  bind_rows(FCS) %>% write_tsv(., paste0(mydir, "/report/", mydir, "_featureCounts.tsv"), col_names = FALSE)
}

# fragment/read-length size (footprint length distribution)
frag_files <- list.files(path = mydirf, pattern = "_fragment_results.tsv")
if (!is_empty(frag_files)) {
  frag_files <- c(paste0(mydirf, "/", frag_files))
  frag <- frag_files %>%
    map(~read_tsv(.x, show_col_types = FALSE)) %>%
    purrr::reduce(full_join, by = "sample")
  sn <- left_join(sn, frag, by = "sample") %>% select(where(~!all(is.na(.))))
}

# riboWaltz P-site/periodicity summary
ribowaltz_file <- list.files(path = mydirf, pattern = "_ribowaltz_psite.tsv", full.names = TRUE)
if (!is_empty(ribowaltz_file)) {
  ribowaltz <- read_tsv(ribowaltz_file[1], show_col_types = FALSE)
  sn <- left_join(sn, ribowaltz, by = "sample") %>% select(where(~!all(is.na(.))))
}

write_tsv(sn, outfile)
