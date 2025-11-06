
library(rmarkdown)
library(docopt)
library(here)

doc <- "Usage: knit_Rmd.R [--help] [--input INPUT] [--output TYPE] [--proj PROJ] [--sample FILE] [--index_map INDEX_MAP] [--seq_date SEQ_DATE] [--index_Sample INDEX_SAMPLE]

Options:
-i --input INPUT        path to rmarkdown
-o --output TYPE        name of output 
-p --proj PROJ          name of project, this is used to name output file
-s --sample FILE        path to sample.yaml file
-m --index_map INDEX_MAP    name of index mapped
-d --seq_date SEQ_DATE      sequence date
-x --index_Sample INDEX_SAMPLE    name of index sample
-h --help               display this help message"

opts <- docopt(doc)

print(opts)

proj         <- opts$proj
res_dir      <- opts$proj
outname      <- opts$output
ymls         <- opts$sample
index_map    <- opts$index_map
index_Sample <- opts$index_Sample
seq_date     <- opts$seq_date
output <- here(paste0(seq_date,  "_", proj,  "_", index_Sample, outname, ".html"))

# Render Rmd
render(
  input       = opts$input,
  output_file = output
)


