
library(rmarkdown)
library(docopt)
library(here)

doc <- "Usage: knit_Rmd.R [--help] [--input INPUT] [--proj PROJ] [--refYaml FILE] [--sample FILE] [--output OUT]

-i --input INPUT    path to rmarkdown
-p --proj PROJ      name of project, this is used to name output file
-o --output OUTPUT  path to directory to write output file
-y --refYaml FILE   path to ref.yaml file
-s --sample FILE    path to sample.yaml file
-h --help           display this help message"

opts <- docopt(doc)

print(opts)

proj    <- opts$proj
res_dir <- opts$output
yml     <- opts$refYaml
ymls    <- opts$sample
output <- here(res_dir, paste0(proj, "_analysis.html"))

# Render Rmd
render(
  input       = opts$input,
  output_file = output
)


