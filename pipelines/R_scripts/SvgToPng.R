#!/usr/bin/env Rscript

# fix text size and postion of svg heatmap files and save as png

library("tidyverse")
library("rsvg")

my_files <- (snakemake@input)
out_files <- (snakemake@output)
my_sizes <- str_split((snakemake@params[["texty"]])," ",simplify = T)

for(i in seq_along(my_files)){
  # Read SVG file
  svg_content <- readLines(my_files[[i]], warn = FALSE)
  svg_content <- paste(svg_content, collapse = "\n")
  my_size <- as.numeric(my_sizes[i])
  
  # Increase font size 
  found_fonts <- str_extract_all(svg_content,"font-size:(\\d+)") %>% unlist() %>% unique() 
  found_fontsize <- found_fonts %>% str_extract_all(.,'(\\d+)') %>% unlist() %>% as.numeric(.) * my_size
  found_fontsize <- floor(found_fontsize)
  modified_svg <- svg_content
  for(j in seq_along(found_fonts)){
    modified_svg <- gsub(
      found_fonts[j],
      paste0("font-size:",found_fontsize[j]),
      modified_svg
    )
  }
  
  # Move text down ~ 1/2 as much as text size
  found_text <- str_extract_all(modified_svg,"<text(.+?)y=\"(\\d+)") %>% unlist() %>% 
    str_replace_all(.,"\\(","\\\\(") %>% 
    str_replace_all(.,"\\)","\\\\)")
  found_texty <- found_text %>% str_split_fixed(.,'y=\"',2) 
  colnames(found_texty) <- c("V1","V2")
  found_texty <- found_texty %>% as_tibble() %>% mutate(V2=as.numeric(V2) + floor(min(found_fontsize)/4)) 
  
  
  modified_svg2 <- modified_svg
  for(j in seq_along(found_text)){
    modified_svg2 <- gsub(
      found_text[j],
      paste0(found_texty[j,1], 'y=\"',found_texty[j,2]),
      modified_svg2
    )
  }
  
  # Create a new temporary SVG file with increased font size and adjusted text position
  temp_svg <- tempfile(fileext = ".svg")
  writeLines(modified_svg2, temp_svg)
  
  # Convert SVG to PNG using rsvg
  rsvg::rsvg_png(temp_svg, file = out_files[[i]], width = 800*my_size, height = 1500)
  
  # Clean up temporary SVG file
  unlink(temp_svg)
}


