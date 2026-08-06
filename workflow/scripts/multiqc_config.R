#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input_tsv   <- strsplit(args[1], " ")[[1]][1] # results.tsv
myproj      <- args[2] # PROJ
seq_date    <- args[3] # SEQ_DATE
index       <- args[4] # INDEXS[0]
output_tsv  <- args[5] # "report.tsv"
output_yaml <- args[6] # "multiqc_config.yaml"

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
})

df <- read_tsv(input_tsv, show_col_types = FALSE)

if ("sample" %in% names(df)) {
  df <- df %>% relocate(sample, .before = everything())
}

headers <- list()

# ---- helpers ----
is_num_percent <- function(x){
  is.character(x) && any(str_detect(x, "%\\)$"))
}

is_percent <- function(x) {
  is.character(x) && any(str_detect(x, "%"))
}

clean_percent <- function(x) {
  as.numeric(str_replace(x, "%", ""))
}

# ---- column processing ----
for (col in names(df)) {
  if (all(is_num_percent(df[[col]]))) {
    df[[col]] <- gsub("\\)$", "", df[[col]]) 
    df <- separate(df,all_of(col),"\\(",into = c(col,paste0(col,"_percent")),convert = T)
  }
}

for (col in names(df)) {
  if (col == "sample") next
  x <- df[[col]]
  
  if (is_percent(x)) {
    df[[col]] <- clean_percent(x)
    headers[[col]] <- list(format = "{:.2f}%", scale = "RdYlGn", min = 0, max = 100)
  } else if (is.numeric(x)) {
    if (all(abs(x - round(x, 0)) < 1e-9, na.rm = TRUE)) {
      headers[[col]] <- list(format = "{:,.0f}", scale = "Blues")
    } else {
      headers[[col]] <- list(format = "{:,.2f}", scale = "Blues")
    }
  } else {
    headers[[col]] <- list(format = "text")
  }
}

# ---- IMAGE SECTION ORDERING & TITLES ----
image_files <- list.files(pattern = paste0("_aligned_",index, "__.*\\.png"), recursive = TRUE, include.dirs = FALSE) %>%
  basename()

prefix <- stringr::str_extract(image_files, "^(543|5L|5|3)")
priority_order <- c("543", "5L", "5", "3")

# Rank each file according to that priority, then reorder image_files by it
ord <- match(prefix, priority_order)
image_files <- image_files[order(ord)]

# Now derive section_ids and clean_titles from the REORDERED image_files
section_ids <- tools::file_path_sans_ext(image_files)
clean_titles <- section_ids %>%
  stringr::str_replace_all("_heatmap$", "") %>%
  stringr::str_replace_all("__.*?bin", "") %>%
  stringr::str_replace_all("^5L", "5") %>%
  stringr::str_replace_all("_", " ") %>%
  stringr::str_trim()

# Establish section priorities inside the custom content module space
section_order_dict <- list()
section_order_dict[["custom_content/Report"]] <- list(order = -1000L)

for (i in seq_along(section_ids)) {
  # Prefixing with custom_content/ forces MultiQC to map the internal ordering correctly
  section_order_dict[[paste0("custom_content/", section_ids[i])]] <- list(order = -1000L + as.integer(i)) 
}

# remove most of the fastqc sub-sections
remove_sections <- list(
  "fastqc-sequence-counts",
  "fastqc-per-base-sequence-quality",
  "fastqc-per-sequence-quality-scores",
  "fastqc-per-base-sequence-content",
  "fastqc-per-sequence-gc-content",
  "fastqc-per-base-n-content",
  "fastqc-sequence-length-distribution",
  "fastqc-sequence-duplication-levels",
  "fastqc-overrepresented-sequences",
  "fastqc-adapter-content"
)

# Loop through and add them to your existing dictionary
for (section in remove_sections) {
  # Convert hyphens to underscores to match MultiQC internal IDs
  clean_id <- gsub("-", "_", section)
  section_order_dict[[clean_id]] <- "remove"
}

header_config <- list()
for (i in seq_along(section_ids)) {
  header_config[[section_ids[i]]] <- list(title = clean_titles[i])
}

# ---- DYNAMICALLY GENERATE SEARCH PATTERNS & CUSTOM DATA ----
sp_list <- list(
  Report = list(fn = basename(output_tsv)),
  bbduk = list(contents = "Executing bbduk", num_lines = 10)
)

custom_data_list <- list(
  Report = list(
    id = "Report",
    section_name = "Overall results",
    description = "Auto-generated summary table",
    plot_type = "table",
    file_format = "tsv",
    headers = headers,
    pconfig = list(id = "Report_table")
  )
)

for (i in seq_along(section_ids)) {
  img_id <- section_ids[i]
  img_file <- image_files[i]
  
  sp_list[[img_id]] <- list(fn = img_file)
  
  custom_data_list[[img_id]] <- list(
    id = img_id,
    section_name = clean_titles[i],
    plot_type = "image"
  )
}

# ---- METHODS TEXT ----
methods_file <- file.path(myproj, "report", paste0(seq_date, "_", myproj, "_pipeline_methods.txt"))

if (file.exists(methods_file)) {
  methods_lines <- readLines(methods_file)
  methods_html  <- paste0(
    "<pre style='font-size:12px; white-space:pre-wrap;'>",
    paste(methods_lines, collapse = "\n"),
    "</pre>"
  )
} else {
  methods_html <- "<p>Methods file not found.</p>"
}

sp_list[["pipeline_methods"]] <- list(fn = basename(methods_file))

custom_data_list[["pipeline_methods"]] <- list(
  id           = "pipeline_methods",
  section_name = "Pipeline Methods",
  description  = "Tool versions and commands used in this pipeline run.",
  plot_type    = "html",
  data         = methods_html
)

section_order_dict[["custom_content/pipeline_methods"]] <- list(order = 1L)


# ---- FULL MULTIQC CONFIG ----
multiqc_config <- list(
  # This wildcard setup forces custom_content to load first, 
  # then catches everything else automatically in its default order.
  top_modules = list(
    "custom_content"
  ),
  
  # Fully qualified sub-names ensure your images sort internally in your exact order
  report_section_order = section_order_dict, 
  
  custom_content_header_config = header_config,
  
  thousandsSep_format = ",",
  decimalPoint_format = ".",
  
  title = myproj,
  subtitle = as.character(seq_date),
  intro_text = "",
  
  skip_versions_section = TRUE,
  skip_generalstats = TRUE,
  ignore_images = FALSE,
  
  sp = sp_list,
  custom_data = custom_data_list
)

# ---- WRITE OUTPUTS ----
yaml_text <- as.yaml(
  multiqc_config,
  handlers = list(logical = verbatim_logical)
)

# Force `title`/`subtitle` to always be emitted as double-quoted YAML
# strings. The yaml package leaves plain-looking scalars unquoted when it
# doesn't think they're ambiguous -- but its notion of "ambiguous" doesn't
# match MultiQC's PyYAML parser, which treats underscore-grouped digits
# (e.g. a run name/date like "2024_01_15" or "240806_A00405") as an
# integer literal when left unquoted. MultiQC then errors expecting a
# string. Explicitly quoting here sidesteps that mismatch regardless of
# what seq_date/myproj look like (plain digits, dashes, underscores, ...).
quote_yaml_line <- function(lines, field, value) {
  escaped <- gsub('(["\\\\])', '\\\\\\1', value)
  ifelse(
    grepl(paste0("^", field, ":"), lines),
    paste0(field, ": \"", escaped, "\""),
    lines
  )
}

lines <- strsplit(yaml_text, "\n", fixed = TRUE)[[1]]
lines <- quote_yaml_line(lines, "title", myproj)
lines <- quote_yaml_line(lines, "subtitle", as.character(seq_date))
yaml_text <- paste(lines, collapse = "\n")

writeLines(yaml_text, output_yaml)

write_tsv(df, output_tsv)
