#!/usr/bin/env Rscript
# fix text size and potions of svg heatmap files and save as png

library("tidyverse")
library("rsvg")
library("xml2")

# ==================== SCRIPT PARAMETERS ====================
# These variables control the behavior of the script and can be easily
# adjusted in one place.

# Text and Title scaling
title_scale <- 1.3         # titles get an extra boost
text_scale <- 1.2          # text gets an extra boost
max_title_len <- 25        # trim long file names

# Positional adjustments
Y_POS_THRESHOLD <- 300     # Threshold for identifying X-axis labels by y-coordinate
LONG_TEXT_THRESHOLD <- 10  # Length threshold to identify titles/long text
Y_SHIFT_DIVISOR <- 1.5     # Divisor for calculating y-axis text shift
X_AXIS_SHIFT <- 5          # Additional y-shift for X-axis labels
GENE_DISTANCE_SHIFT <- 5   # Additional y-shift for "gene distance" text

# Line and Tick scaling
TICK_MULTIPLIER <- 4       # Multiplier for tick length
MAX_STROKE_WIDTH_STYLE <- 500 # Maximum stroke width for styled lines
CAP_WIDTH_STYLE <- 5         # Cap value for styled lines
MAX_STROKE_WIDTH_ATTR <- 5 # Cap value for attributed lines
CAP_WIDTH_ATTR <- 5        # Cap value for attributed lines

# Padding and output
PADDING_MULTIPLIER <- 20   # Multiplier for padding calculation
PNG_WIDTH <- 800           # Base width for the output PNG
PNG_HEIGHT_BASE <- 1500    # Base height for the output PNG
PNG_HEIGHT_MULTIPLIER <- 300 # Multiplier for calculating additional PNG height


# ==================== MAIN SCRIPT LOGIC ====================

args = commandArgs(trailingOnly=TRUE)
my_files <- strsplit(args[1], " ")[[1]]
out_files <- strsplit(args[3], " ")[[1]]
my_sizes <- strsplit(args[2], " ")[[1]]


for(i in seq_along(my_files)){
  # Added tryCatch block to gracefully handle potential XML parsing errors
  tryCatch({
    # Load SVG
    doc <- read_xml(my_files[[i]])
    my_size <- as.numeric(my_sizes[i])
    tick_length <- TICK_MULTIPLIER * my_size  # new tick length
    
    ## ---- FIND PLOT BOUNDS FOR CENTERING ----
    # Find the horizontal center of the main plot area by looking at the rect elements.
    rects <- xml_find_all(doc, ".//rect")
    if (length(rects) > 0) {
      x_coords <- as.numeric(xml_attr(rects, "x"))
      # Filter out potential background rectangles by taking a quantile
      x_coords <- x_coords[x_coords > 0] 
      min_x <- min(x_coords, na.rm = TRUE)
      max_x <- max(x_coords, na.rm = TRUE)
      plot_center_x <- (min_x + max_x) / 2
    } else {
      plot_center_x <- NA
    }
    
    ## ---- TEXT ----
    # Scale all text sizes
    doc_clean <- xml_ns_strip(doc)
    texts <- xml_find_all(doc_clean, ".//text")
    
    # Define X position threshold for Y-axis labels (you may need to adjust this)
    
    for (t in texts) {
      style <- xml_attr(t, "style")
      content <- xml_text(t)
      
      # Check if the text is an X-axis label by its y-coordinate
      is_x_axis_label <- FALSE
      y_pos <- as.numeric(xml_attr(t, "y"))
      if (!is.na(y_pos) && y_pos > Y_POS_THRESHOLD) {
        is_x_axis_label <- TRUE
      }
      
      # Check if the text is a Y-axis label by its text-anchor
      is_y_axis_label <- FALSE
      text_anchor <- xml_attr(t, "text-anchor")
      if (!is.na(text_anchor) && text_anchor == "end") {
        is_y_axis_label <- TRUE
      }
      
      # Extract font-size from style string
      size_match <- str_match(style, "font-size:\\s*([0-9.]+)px")
      
      if (!is.na(size_match[2])) {
        size_num <- as.numeric(size_match[2])
        
        # Determine scaling factor
        is_title <- (nchar(content) > LONG_TEXT_THRESHOLD & !grepl("^-?\\d*\\.?\\d+$", content) & !str_detect(content, "gene distance"))
        
        if (is_title) {
          # This is a title or long text
          style <- paste(style, "font-weight: bold;", sep = " ")
          # Explicitly ensure proper text anchoring
          xml_set_attr(t, "text-anchor", "middle")
          
          newsize <- size_num * my_size * title_scale
          if (nchar(content) > max_title_len) {
            # First remove whitespace, dashes, and underscores
            cleaned_content <- gsub("[\\s_-]", "", content)
            if (nchar(cleaned_content) > max_title_len) {
              cleaned_content <- paste0(substr(cleaned_content, 1, max_title_len-3), "…")
            }
            xml_text(t) <- cleaned_content
          }
          
          # Update the font size in the style
          new_style <- sub("font-size:\\s*[0-9.]+px", 
                           paste0("font-size: ", newsize, "px"), 
                           style)
          xml_set_attr(t, "style", new_style)
          
        } else {
          # This is a regular label
          newsize <- size_num * my_size * text_scale
          
          # Calculate the Y-axis shift
          y_shift <- (newsize - size_num) / Y_SHIFT_DIVISOR
          
          # Get the original Y coordinate
          y_orig <- as.numeric(xml_attr(t, "y"))
          
          # Check if this is a cluster label that needs rotation adjustment
          is_cluster_label <- str_detect(content, "cluster_")
          
          if (!is.na(y_orig)) {
            # Apply different shifts based on label type
            if (is_y_axis_label) {
              # No shift for Y-axis labels (they have text-anchor="end")
              new_y <- y_orig
            } else if (is_x_axis_label) {
              # Apply full shift for X-axis labels
              new_y <- y_orig + y_shift + X_AXIS_SHIFT
            } else {
              # Apply normal shift for other labels
              new_y <- y_orig + y_shift
            }
            
            # Handle cluster labels rotation
            if (is_cluster_label) {
              # Get current transform attribute
              current_transform <- xml_attr(t, "transform")
              x_pos <- as.numeric(xml_attr(t, "x"))
              
              if (!is.na(x_pos)) {
                # Change rotation from -90 to -45 degrees (or 0 for horizontal)
                new_transform <- paste0("rotate(-45 ", x_pos, " ", new_y, ")")
                # For horizontal: new_transform <- paste0("rotate(0 ", x_pos, " ", new_y, ")")
                
                xml_set_attr(t, "transform", new_transform)
                
                # Adjust positioning to prevent overlap - move left slightly
                xml_set_attr(t, "x", as.character(x_pos - 40))  # Adjust as needed
                PADDING_MULTIPLIER <- 40
              }
            }
            
            # Apply extra shift for "gene distance"
            if (str_detect(content, "gene distance")) {
              new_y <- new_y + GENE_DISTANCE_SHIFT
            }
            
            xml_set_attr(t, "y", as.character(new_y))
          }
          
          # Replace font-size in style string
          new_style <- sub("font-size:\\s*[0-9.]+px", 
                           paste0("font-size: ", newsize, "px"), 
                           style)
          xml_set_attr(t, "style", new_style)
        }
      }
    }
    
    ## ---- LINES ----
    # Increase line size
    # Convert xml2 object back to string 
    modified_svg2 <- as.character(doc_clean)
    
    # Find and scale stroke-width in style attributes (for tick marks)
    style_stroke_pattern <- 'stroke-width:\\s*(\\d+(?:\\.\\d+)?)'
    stroke_matches <- str_extract_all(modified_svg2, style_stroke_pattern)[[1]]
    
    if(length(stroke_matches) > 0) {
      stroke_matches <- unique(stroke_matches)
      
      for(stroke_match in stroke_matches) {
        # Extract current width value
        current_width <- as.numeric(str_extract(stroke_match, "\\d+(?:\\.\\d+)?"))
        
        # Scale stroke width for tick marks (avoid making them too thick)
        new_width <- current_width * my_size
        if(new_width > MAX_STROKE_WIDTH_STYLE) new_width <- CAP_WIDTH_STYLE
        
        # Replace in SVG
        new_stroke <- paste0("stroke-width: ", new_width)
        modified_svg2 <- gsub(stroke_match, new_stroke, modified_svg2, fixed = TRUE)
      }
    }
    
    # Also handle any remaining attribute-style stroke-width (backup)
    attr_stroke_pattern <- 'stroke-width="([^"]*)"'
    attr_matches <- str_match_all(modified_svg2, attr_stroke_pattern)[[1]]
    
    if(nrow(attr_matches) > 0) {
      for(m in 1:nrow(attr_matches)) {
        full_attr <- attr_matches[m, 1]
        current_width <- as.numeric(attr_matches[m, 2])
        
        new_width <- current_width * my_size
        if(new_width > MAX_STROKE_WIDTH_ATTR) new_width <- CAP_WIDTH_ATTR
        
        new_attr <- paste0('stroke-width="', new_width, '"')
        modified_svg2 <- gsub(full_attr, new_attr, modified_svg2, fixed = TRUE)
      }
    }
    # Convert back to xml2 object
    doc_clean <- read_xml(modified_svg2)
    
    ## ---- TICKS ----
    ticks <- xml_find_all(doc_clean, ".//path[@id]")
    for (tk in ticks) {
      d <- xml_attr(tk, "d")
      if (!is.na(d)) {
        d <- sub("M 0 0  L 0 [0-9.]+", paste0("M 0 0  L 0 ", tick_length), d)
        d <- sub("M 0 0  L -[0-9.]+ 0", paste0("M 0 0  L -", tick_length, " 0"), d)
        xml_set_attr(tk, "d", d)
      }
    }
    
    ## ---- padding ----
    # Convert xml2 object back to string to use your working method
    modified_svg2 <- as.character(doc_clean)
    
    # Add padding to SVG by modifying viewBox and adding background
    padding_px <- PADDING_MULTIPLIER * my_size
    
    # Extract current viewBox or create one
    viewbox_match <- str_match(modified_svg2, 'viewBox="([^"]*)"')
    if(!is.na(viewbox_match[1,1])) {
      # Parse existing viewBox
      viewbox_values <- as.numeric(strsplit(viewbox_match[1,2], "\\s+")[[1]])
      orig_x <- viewbox_values[1]
      orig_y <- viewbox_values[2] 
      orig_width <- viewbox_values[3]
      orig_height <- viewbox_values[4]
    } else {
      # Extract width/height attributes if no viewBox
      width_match <- str_extract(modified_svg2, 'width="([^"]*)"')
      height_match <- str_extract(modified_svg2, 'height="([^"]*)"')
      orig_x <- 0
      orig_y <- 0
      orig_width <- as.numeric(str_extract(width_match, "\\d+"))
      orig_height <- as.numeric(str_extract(height_match, "\\d+"))
    }
    
    # Calculate new dimensions with padding
    new_x <- orig_x - padding_px
    new_y <- orig_y - padding_px
    new_width <- orig_width + (2 * padding_px)
    new_height <- orig_height + (2 * padding_px)
    new_viewbox <- paste(new_x, new_y, new_width, new_height)
    
    # Update or add viewBox
    if(!is.na(viewbox_match[1,1])) {
      modified_svg3 <- gsub('viewBox="[^"]*"', paste0('viewBox="', new_viewbox, '"'), modified_svg2)
    } else {
      modified_svg3 <- gsub('<svg', paste0('<svg viewBox="', new_viewbox, '"'), modified_svg2)
    }
    
    # Add white background rectangle
    bg_rect <- paste0('<rect x="', new_x, '" y="', new_y, '" width="', new_width, 
                      '" height="', new_height, '" fill="white"/>')
    
    # Insert background rectangle right after opening svg tag
    modified_svg3 <- gsub('(<svg[^>]*>)', paste0('\\1\n', bg_rect), modified_svg3)
    
    # Convert back to xml2 object
    doc_clean <- read_xml(modified_svg3)
    
    
    ## ---- SAVE & CONVERT ----
    tmp_svg <- tempfile(fileext = ".svg")
    write_xml(doc_clean, tmp_svg)
    rsvg::rsvg_png(tmp_svg, file = out_files[[i]], 
                   width = PNG_WIDTH * my_size, 
                   height = PNG_HEIGHT_BASE + (PNG_HEIGHT_MULTIPLIER * my_size))
    
    # Clean up temp file
    unlink(tmp_svg)
    
  }, error = function(e) {
    message("Error processing file ", my_files[[i]], ": ", e$message)
  })
  
}
