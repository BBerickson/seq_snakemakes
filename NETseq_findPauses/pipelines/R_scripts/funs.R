

# functions

#' Plot color palette
#' 
#' @param cols_in Character vector containing colors to plot
#' @param ... Additional arguments to pass to ggplot
#' @return ggplot object
#' @export
plot_color_palette <- function(cols_in, ...) {
  
  if (is.null(names(cols_in))) {
    cols_in <- set_names(cols_in)
  }
  
  col_df <- tibble::tibble(color = factor(names(cols_in), names(cols_in)))
  
  res <- ggplot2::ggplot(col_df, aes(x = "color", fill = color)) +
    ggplot2::geom_bar(...) +
    ggplot2::scale_fill_manual(values = cols_in) +
    ggplot2::theme_void()
  
  res
}

#' Format p-values for plotting
#' @noRd
.format_pvalue <- function(p, digits = 1, cutoffs = NULL) {
  
  # Set p label based on vector of cutoffs
  if (!is.finite(p)) return(as.character(NA))
  
  if (!is.null(cutoffs)) {
    if (any(duplicated(cutoffs))) {
      cli::cli_abort("Cutoff values for p_label must be unique.")
    }
    
    # Set default labels when not provided by user
    if (is.null(names(cutoffs))) {
      cutoffs <- sort(cutoffs, decreasing = TRUE)
      
      names(cutoffs) <- purrr::imap_chr(
        cutoffs, ~ paste0(rep("*", .y), collapse = "")
      )
    }
    
    cutoffs <- sort(cutoffs)
    p_label <- as.character(NA)
    
    for (val in names(cutoffs)) {
      if (p < cutoffs[val]) {
        p_label <- val
        
        break()
      }
    }
    
    # Treat "value" as a keyword that will allow user to display actual
    # p-value for a certain cutoff
    # All custom labels need to be wrapped in quotes for parsing
    if (!identical(p_label, "value")) {
      if (!is.na(p_label)) p_label <- paste0("\'", p_label, "\'")
      
      return(p_label)
    }
  }
  
  # Format p-value label
  if (p >= 0.1) return(as.character(round(p, 1)))
  
  p <- scales::label_scientific(digits = digits)(p)
  
  ex <- str_extract_all(p, "[+\\-][0-9]+$")
  
  p <- sub(paste0("\\", ex, "$"), "", p)
  
  ex <- as.numeric(ex)
  ex <- as.character(ex)
  
  p <- sub("e", "*x*10^", p)
  p <- paste0(p, ex)
  
  p
}

#' Create labeller function to add cell n labels
#' 
#' @param sobj_in Seurat object or data.frame containing plot data.
#' @param lab_col meta.data column containing cell groups.
#' @param nm Should cell group be included in label.
#' @param sep Separator to use for labels.
#' @return Labeller function
#' @export
get_nlab_fun <- function(df_in, lab_col = NULL, nm = TRUE, sep = "\n") {
  
  if (!is.null(lab_col)) {
    df_in <- df_in %>%
      group_by(!!sym(lab_col))
  }
  
  labs <- df_in %>%
    summarize(n = n(), .groups = "drop") %>%
    distinct() %>%
    arrange(desc(n)) %>%
    mutate(
      n = format(n, big.mark = ",", scientific = FALSE),
      n = str_trim(n),
      n = str_c("n = ", n)
    )
  
  if (nm && !is.null(lab_col)) {
    labs <- labs %>%
      mutate(n = str_c(!!sym(lab_col), sep, "(", n, ")"))
  }
  
  if (is.null(lab_col)) {
    return(pull(labs, "n"))
  }
  
  labs <- set_names(
    x  = pull(labs, "n"),
    nm = pull(labs, lab_col)
  )
  
  # res <- as_labeller(labs)
  
  res <- function(x) labs[x]
  
  res
}

# Remove paths and fastq info from sample names
trim_path <- function(input_df, rm_str = NULL, file_col = "sample") {
  
  file_sym <- sym(file_col)
  
  res <- input_df %>%
    mutate(!!file_sym := base::basename(!!file_sym))
  
  if (!is.null(rm_str)) {
    for (i in rm_str) {
      res <- res %>%
        mutate(!!file_sym := str_remove(!!file_sym, i))
    }
  }
  
  res
}

# Function to remove string common to all sample names  
remove_com_str <- function(input_str, input_names = NULL, sep = "_") {
  
  if (is.null(input_names)) {
    input_names <- input_str
  }
  
  input_names <- input_names %>%
    unique()
  
  # Unique strings in keys
  uniq_strs <- input_names %>%
    str_split(sep) %>%
    unlist() %>% 
    unique()
  
  uniq_strs <- uniq_strs[uniq_strs != ""]
  
  # Identify strings that are present in every key
  common_strs <- NULL
  
  for (i in uniq_strs) {
    str_num <- input_names %>%
      str_detect(i) %>%
      (function(x) length(x[x == T]))
    
    if (str_num == length(input_names)) {
      common_strs <- c(common_strs, i)
    }
  }
  
  if (is.null(common_strs)) {
    return(input_str)
  }
  
  # Duplicate sep that needs to be removed
  internal_sep <- str_c(sep, sep)
  
  common_strs <- common_strs %>%
    purrr::reduce(str_c, sep = "|")
  
  res <- input_str %>%
    str_remove(common_strs) %>%
    str_remove_all(internal_sep) %>%
    str_remove(str_c("^", sep)) %>%
    str_remove(str_c(sep, "$"))
  
  res
}

# Remove strings that are common in all sample names
remove_com_key <- function(input_df, key_col = "sample", sep = "_") {
  
  # Get unique keys
  input_names <- input_df[, key_col] %>%
    unique() %>%
    unlist() %>%
    as.character()
  
  # Remove common string from keys
  key_sym <- sym(key_col)
  
  res <- input_df %>%
    mutate(!!key_sym := remove_com_str(!!key_sym, input_names, sep = sep))
  
  res
}

# Function to retrieve unique gene names
get_uniq_names <- function(input_df) {
  input_df %>%
    dplyr::select(name) %>%
    unique()
}

# Shorten sample names
shorten_names <- function(nm, ...) {
  if (!is_empty(list(...))) {
    nm <- nm %>%
      str_remove_all(...)
  }
  
  nm %>%
    str_remove_all("^NF_|_mNETseq") %>%
    str_remove_all("[\\-_]+$")
}

# Set axis breaks as the limits
breaks_limits <- function(x) {
  
  round_down <- function(x, to = 10) {
    # res <- to * (res %/% to + as.logical(res %% to))  # this rounds up
    res <- abs(x)
    res <- to * (res %/% to)
    
    res[x < 0] <- res[x < 0] * -1
    
    res
  }
  
  count_digits <- function(x) {
    res <- floor(log10(abs(x))) + 1
    res[!is.finite(res)] <- 0
    res
  }
  
  # Remove padding from limits
  pad_f  <- 0.05
  no_pad <- diff(x) / (pad_f * 2 + 1)
  pad    <- no_pad * pad_f
  
  x[1] <- x[1] + pad
  x[2] <- x[2] - pad
  
  # Round down limits
  to         <- count_digits(x)
  to         <- to - 2
  to[to < 0] <- 0
  to         <- 10 ^ to
  
  res <- round_down(x, to)
  
  res[x < 10] <- round(x[x < 10], 1)
  
  res
}

# Extract distance genes are separated using gene list file name
extract_sep_info <- function(x) {
  res <- x %>%
    basename() %>%
    str_extract("[0-9]+(?=ksep)")
  
  res
}

# Format sample names for plotting
format_sample_names <- function(df_in, key_vec = sam_lnms, add_grp = FALSE) {
  if (add_grp) {
    df_in <- df_in %>%
      separate(
        sample,
        into   = c(NA, "group"),
        sep    = "_(?=[^_]*$)",
        extra  = "merge",
        remove = FALSE
      )
  }
  
  res <- df_in %>%
    mutate(sample = key_vec[sample]) %>%
    separate(sample, into = c("treat", "rep"), sep = "_(?=[^_]*$)", remove = FALSE)
  
  res
}

# Create QC bar graphs
create_qc_bars <- function(df_in, grp_df, grp_lvls = NULL, sam_lvls = NULL, met_lvls = NULL, plot_cols = NULL,
                           lab_metric = NULL, lab_div = 1e6, lab_unit = "million", lab_size = 9, n_rows = 1, ...) {
  
  # Calculate metrics
  clmns <- "sample"
  
  if ("grp" %in% names(df_in)) {
    clmns <- c("grp", clmns)
  }
  
  res <- df_in %>%
    group_by(!!!syms(clmns)) %>%
    mutate(
      reads = sum(value),
      frac  = value / reads
    ) %>%
    ungroup()
  
  # Create bar label
  if (!is.null(lab_metric)) {
    if (any(!lab_metric %in% res$metric)) {
      stop(str_c("Metrics (", str_c(lab_metric, collapse = ", "), ") not found in data.frame."))
    }
    
    lab_df <- res %>%
      dplyr::filter(metric %in% lab_metric) %>%
      group_by(!!!syms(clmns)) %>%
      summarize(
        lab     = sum(value) / lab_div,
        lab     = round(lab, 1),
        lab     = str_c(lab, " ", lab_unit),
        .groups = "drop"
      ) %>%
      ungroup() %>%
      dplyr::select(all_of(clmns), lab)
    
    res <- res %>%
      left_join(lab_df, by = clmns)
  }
  
  # Add group names
  if (!"grp" %in% names(res) && !is.null(grp_df)) {
    if (!all(unique(res$sample) %in% grp_df$sample)) {
      warning("Not all samples present in data are included in grp_df, some samples not shown.")
    }
    
    res <- res %>%
      left_join(grp_df, by = "sample") %>%
      dplyr::filter(!is.na(grp))
  }
  
  # Set factor levels
  if (!is.null(grp_lvls)) {
    if (!all(unique(res$grp) %in% grp_lvls)) {
      warning("Not all groups present in data are included in grp_lvls, some groups not shown.")
    }
    
    res <- res %>%
      dplyr::filter(grp %in% grp_lvls) %>%
      mutate(grp = fct_relevel(grp, grp_lvls))
  }
  
  if (!is.null(sam_lvls)) {
    if (!all(unique(res$sample) %in% names(sam_lvls))) {
      warning("Not all samples present in data are included in sam_lvls.")
    }
    
    res <- res %>%
      dplyr::filter(sample %in% names(sam_lvls)) %>%
      mutate(
        sample = recode(sample, !!!sam_lvls),
        sample = fct_relevel(sample, sam_lvls)
      )
  }
  
  if (!is.null(met_lvls)) {
    if (!all(unique(res$metric) %in% met_lvls)) {
      warning("Not all metrics present in data are included in met_lvls, some metrics not shown.")
    }
    
    res <- res %>%
      dplyr::filter(metric %in% met_lvls) %>%
      mutate(metric = fct_relevel(metric, met_lvls))
  }
  
  # Create bar graphs
  res <- res %>%
    ggplot(aes(sample, frac, fill = metric)) +
    geom_col(color = "white", size = 0.3, ...) +
    
    guides(fill = guide_legend(nrow = 3)) +
    
    theme_info +
    theme(
      plot.margin     = margin(5, 5, 5, 40),
      legend.position = "top",
      legend.title    = element_blank(),
      strip.text      = element_text(size = 10),
      axis.title      = element_blank(),
      axis.text.x     = element_text(angle = 45, hjust = 1)
    )
  
  # Split into facets
  if ("grp" %in% names(df_in) || !is.null(grp_df)) {
    res <- res +
      facet_wrap(~ grp, scales = "free_x", nrow = n_rows)
  }
  
  # Add plot labels
  if (!is.null(lab_metric)) {
    res <- res +
      geom_text(
        aes(y = 0.1, label = lab),
        check_overlap = TRUE,
        angle         = 90,
        hjust         = 0,
        size          = lab_size / .pt
      )
  }
  
  # Set plot colors
  if (!is.null(plot_cols)) {
    res <- res +
      scale_fill_manual(values = plot_cols)
  }
  
  res
}

#' Helper to remove genes that are not shared between categorical variables
#'
#' @param df_in tibble or data.frame
#' @param clmns Columns to use for filtering genes. Genes will only be kept if
#' data is present for all variables stored in each column.
#' @param gene_clmn Column containing gene names.
#' @return tibble
#' @export
get_overlapping_genes <- function(df_in, clmns, gene_clmn = "name") {
  
  df <- df_in %>% 
    distinct(!!!syms(c(gene_clmn, clmns)))
  
  clmns %>%
    walk(~ {
      uniq_grps <- unique(df[[.x]])
      grp_clmns <- clmns[clmns != .x]
      
      df <<- df %>%
        group_by(!!!syms(c(gene_clmn, grp_clmns))) %>%
        dplyr::filter(all(uniq_grps %in% !!sym(.x))) %>%
        ungroup()
    })
  
  if (nrow(df) == 0) {
    stop("There are no genes remaining after filtering")
  }
  
  df <- df %>%
    distinct(!!sym(gene_clmn))
  
  res <- df_in %>%
    semi_join(df, by = gene_clmn)
  
  res <- res %>%
    check_group_size(clmns, warn = TRUE, quiet = TRUE)
  
  res
}

#' Helper to check that categorical variables appear an equal number of times
#'
#' @param df_in tibble or data.frame
#' @param clmns Columns to check. An error will be thrown if variables stored
#' in each column are not present an equal number of times. If set to NULL all
#' columns of type character are checked.
#' @param warn If set to TRUE, instread of throwing an error, a warning will be
#' displayed.
#' @param quiet If set to TRUE, messages and warnings will not be shown.
#' @return Messages, warnings, and errors describing the results
#' @export
check_group_size <- function(df_in, clmns = NULL, warn = FALSE, quiet = FALSE) {
  
  .check_clmns <- function(df, cls) {
    res <- cls %>%
      map(~ {
        lens <- unique(table(df[[.x]]))
        
        if (length(lens) != 1) {
          .x
          
        } else if (!quiet) {
          message(.x, ": ", lens)
        }
      })
    
    res <- purrr::discard(res, is.null)
  }
  
  # Identify character columns
  other_clmns <- df_in %>%
    select_if(is.character) %>%
    names()
  
  other_clmns <- other_clmns[!other_clmns %in% clmns]
  
  # Identify other columns to warn
  if (!quiet && !is.null(clmns)) {
    warn <- df_in %>%
      .check_clmns(other_clmns)
    
    if (length(warn) > 0) {
      warn <- str_c(warn, collapse = ", ")
      
      warning("Groups in other columns (", warn, ") of type character are not present an equal number of times")
    }
  }
  
  # Identify bad columns
  if (is.null(clmns)) {
    clmns <- other_clmns
  }
  
  bad <- df_in %>%
    .check_clmns(clmns)
  
  if (length(bad) > 0) {
    bad <- str_c(bad, collapse = ", ")
    
    err_msg <- str_c("Groups in ", bad, " are not present an equal number of times")
    
    if (!warn) {
      stop(err_msg)
    }
    
    warning(err_msg)
  }
  
  df_in
}

#' Helper to load bed files
#'
#' @param path Path to bed file to load
#' @param col_names Column names to use for bed file
#' @param genes Gene list to use for filtering bed file
#' @param join_by Columns to use for joining with genes
#' @param ... Additional sample labels to add as new columns in the tibble
#' @return tibble
#' @export
load_bed <- function(path, col_names, genes = NULL, join_by = "name", ...) {
  
  res <- path %>%
    vroom(delim = "\t", col_names = col_names)
  
  if (!is.null(genes) && !is.null(join_by)) {
    res <- res %>%
      semi_join(genes, by = join_by)
  }
  
  res <- res %>%
    mutate(...)
  
  res
}

#' Helper to load and merge bed files
#'
#' @param prfxs Vector containing file prefixes. The file prefixes should
#' also correspond to the sample directory. The combination of the file
#' prefixes and suffixes must produce a complete file name.
#' @param sfxs Named vector containing file suffixes. Names should correspond
#' to unique sample type, which will get added as a new column in the tibble.
#' The combination of the file prefixes and suffixes must produce a complete
#' file name.
#' @param paths Vector of directory paths where files can be found. If a single
#' path is provided, this will be used for all files.
#' @param group Group name for the file. Group name will get added as a new
#' column in the tibble.
#' @param genes Gene list to use for filtering bed file.
#' @param file_out Path to output file to save tibble.
#' @param overwrite Overwrite file_out if it exists.
#' @param col_names Column names to use for bed file.
#' @param expr_quartiles Should an additional column be added containing gene
#' expression quartiles.
#' @param ... Additional arguments to pass to merge_wins.
#' @return tibble
#' @export
load_merge_wins <- function(prfxs, sfxs, paths, group, genes = NULL, file_out = NULL,
                            overwrite = FALSE, col_names, expr_quartiles = TRUE, ...) {
  
  if (!is.null(file_out) && file.exists(file_out) && !overwrite) {
    return(vroom(file_out))
  }
  
  if (length(paths) != length(prfxs) && length(paths) != 1) {
    stop("paths must be either length 1 or the same length as prfxs")
  }
  
  if (length(paths) == 1) {
    paths <- rep(paths, length(prfxs))
  }
  
  paths <- set_names(paths, prfxs)
  
  paths <- crossing(
    group  = group,
    sample = prfxs,
    type   = names(sfxs)
  ) %>%
    mutate(
      path = str_c(sample, sfxs[type]),
      path = file.path(paths[sample], path)
    )
  
  # Load bed files
  wins <- paths %>%
    pmap_dfr(~ {
      args <- list(...)
      
      args$path %>%
        load_bed(
          col_names = col_names,
          genes     = genes
        ) %>%
        mutate(
          group  = args$group,
          sample = args$sample,
          type   = args$type,
          file   = basename(args$path)
        )
    })
  
  # Merge window bed files so all files include the same genes
  res <- wins %>%
    merge_wins(..., groups = c("sample", "group"))
  
  # Divide genes based on expression
  if (expr_quartiles) {
    qrts <- res %>%
      group_by(sample, group, name, type) %>%
      summarize(qrt = sum(counts), .groups = "drop") %>%
      group_by(group, name, type) %>%                     # Average for all samples
      summarize(qrt = mean(qrt), .groups = "drop") %>%
      group_by(type) %>%
      mutate(qrt = ntile(qrt, 4)) %>%
      ungroup()
    
    res <- res %>%
      left_join(qrts, by = c("group", "name", "type"))
  }
  
  # Check for NAs
  if (!all(complete.cases(res))) {
    stop("NAs present in merged data.frame.")
  }
  
  # Check for duplicate rows
  if (nrow(res) != nrow(distinct(res))) {
    stop("Duplicate rows present in merged data.frame.")
  }
  
  # Write data.frame
  if (!is.null(file_out)) {
    vroom_write(res, file_out, delim = "\t")
  }
  
  res
}

# Helper to load gene region bed files
load_region_beds <- function(sam, reg, genes = NULL, dat_dir = res_dir, prfx = "_",
                             col_names = c("chrom", "start", "end", "name", "score", "strand", "counts")) {
  
  create_path <- function(bed_dir, type, sfx) {
    str_c(dat_dir,"/", bed_dir,"/", str_c(sam, type, reg, sfx))
  }
  
  files <- c(
    NET       = create_path("metaplot_beds", "_", "_S.bed.gz"),
    `p-reads` = create_path("beds/pauses", str_c(pause_prfx, "pause_reads_"), ".bed.gz"),
    pauses    = create_path("beds/pauses", str_c(pause_prfx, "pauses_"), ".bed.gz")
  )
  
  res <- files %>%
    imap_dfr(~ {
      df <- .x %>%
        vroom(delim = "\t", col_names = col_names) %>%
        mutate(
          sample = sam,
          region = reg,
          type   = .y
        )
      
      if (!is.null(genes)) {
        df <- df %>%
          semi_join(genes, by = "name")
      }
      
      df
    })
  
  res
}

#' Calculate mean signal
#'
#' @param df_in data.frame containing signal for each gene.
#' @param rel Calculate mean relative signal. This expresses the signal
#' relative to total signal in the gene (window signal / sum(total gene signal))
#' @param win_col Column containing window IDs.
#' @param flip_type Express signal as negative values for certain samples. This
#' should be a named vector indicating the value to look for when determining
#' which rows to multiply by -1. The vector name indicates the column.
calc_mean_signal <- function(df_in, rel = FALSE, win_col = "win_dist", flip_type = c("type" = "anti-sense"),
                             grps = c("sample", "group", "type")) {
  
  res <- df_in
  
  # Calculate relative signal
  if (rel) {
    rel_grps <- grps
    
    if (!is.null(flip_type)) {
      rel_grps <- grps[grps != names(flip_type)]
    }
    
    res <- res %>%
      group_by(!!!syms(rel_grps), name) %>%
      mutate(counts = counts / sum(counts)) %>%
      ungroup()
  }
  
  res <- res %>%
    group_by(!!!syms(c(grps, win_col))) %>%
    summarize(
      counts  = mean(counts),
      n       = n_distinct(name),
      n       = str_c("n = ", comma(n)),
      .groups = "drop"
    ) %>%
    mutate(counts = counts * 1000)
  
  # Flip negative strand values
  if (!is.null(flip_type)) {
    res <- res %>%
      mutate(counts = ifelse(!!sym(names(flip_type)) == unname(flip_type), counts * -1, counts))
  }
  
  res 
}

# Helper to create meta plots
create_meta <- function(df_in, x, y, color = NULL, alpha = NULL, plot_clrs = NULL, plot_alphs = c(0.5, 1), plot_lvls = NULL,
                        size = 1, h_line = 0, v_line = 0, vline_clr = "black", vline_type = 2, n_lab = NULL,
                        n_lab_pos = c(Inf, Inf), se_clmn = NULL, ...) {
  
  # Check data.frame for columns
  clmns <- c(x, y, color, alpha, n_lab, se_clmn)
  
  need <- clmns[!clmns %in% colnames(df_in)]
  
  if (length(need) > 0) {
    stop("Some columns are not present in the provided data.frame: ", str_c(need, collapse = ", "))
  }
  
  # Calculate standard error
  # By providing se_clmn, this allows a check on the expected number of values
  # used for the calculation
  # The number of values used to calculate SE should equal the number of unique
  # values in se_clmn
  res <- df_in
  
  if (!is.null(color)) {
    res <- res %>%
      mutate(!!sym(color) := fct_relevel(!!sym(color), unique(plot_lvls)))
  }
  
  if (!is.null(se_clmn)) {
    res <- res %>%
      group_by(!!!syms(c(x, color, alpha, n_lab)))
    
    res <- res %>%
      summarize(
        se_min    = min(!!sym(y)),
        se_max    = max(!!sym(y)),
        !!sym(y) := mean(!!sym(y)),
        .groups   = "drop"
      )
  }
  
  # Create base plot
  gg_args <- list(
    x = x,
    y = y
  )
  
  rib_args <- list(
    x    = x,
    ymin = "se_min",
    ymax = "se_max"
  )
  
  if (!is.null(color)) {
    gg_args$color <- color
    rib_args$fill <- color
  }
  
  if (!is.null(alpha)) {
    gg_args$alpha     <- alpha
    rib_args$linetype <- alpha
  }
  
  res <- res %>%
    ggplot(do.call(aes, syms(gg_args)))
  
  if (!is.null(se_clmn)) {
    res <- res +
      geom_ribbon(
        do.call(aes, syms(rib_args)),
        show.legend = FALSE,
        alpha       = 0.25,
        color       = NA
      )
  }
  
  res <- res +
    geom_vline(xintercept = v_line, color = vline_clr, linetype = vline_type) +
    geom_hline(yintercept = h_line, linewidth = 0.25) +
    geom_line(linewidth = size, key_glyph = draw_key_point, ...)
  
  # Add n labels
  # Check that there is only one n label present
  if (!is.null(n_lab)) {
    n_n_lab <- df_in %>%
      pull(n_lab) %>%
      n_distinct()
    
    if (n_n_lab > 1) {
      stop("There are multiple n labels present in ", n_lab)
    }
    
    res <- res +
      geom_text(
        mapping       = aes(x = n_lab_pos[1], y = n_lab_pos[2], label = !!sym(n_lab)),
        color         = "black",
        size          = txt_pt / .pt,
        alpha         = 1,
        check_overlap = TRUE,
        hjust         = 1.4,
        vjust         = 2
      )
  }
  
  # Adjust plot aesthetics
  res <- res +
    guides(
      color = guide_legend(override.aes = list(size = 3), ncol = 1),
      alpha = "none"
    ) +
    
    theme_info +
    theme(
      plot.title      = element_text(face = "plain"),
      legend.position = c(0.01, 0.9),
      legend.title    = element_blank(),
      legend.text     = element_text(size = ttl_pt),
      axis.title.x    = element_blank(),
    )
  
  if (!is.null(plot_clrs)) {
    res <- res +
      scale_color_manual(values = plot_clrs) +
      scale_fill_manual(values = plot_clrs)
  }
  
  if (!is.null(plot_alphs)) {
    res <- res +
      scale_alpha_manual(values = plot_alphs)
  }
  
  res
}

# Helper to add x-axis breaks
add_breaks <- function(breaks = seq(-10, 10), zero_lab = "0 kb") {
  x_fun <- function(x) {
    if_else(round(x, 10) == 0, zero_lab, str_c(x, " kb"))  # x-axis zero is actually stored as a really small number >0
  }
  
  scale_x_continuous(
    labels = x_fun,
    breaks = breaks,
  )
}

# Helper to create 5' and 3' meta plot panels
create_meta_fig <- function(df_5, df_3, color, ylim_3, sams = NULL, grp = NULL,
                            plot_clrs, se_clmn = NULL, file = NULL, dev = "png",
                            dims = c(10, 6), plot_ttl = "mean signal (RPKM)",
                            return_list = FALSE, ...) {
  
  # Filter samples
  if (!is.null(sams)) {
    df_5 <- df_5 %>%
      dplyr::filter(!!sym(color) %in% sams)
    
    df_3 <- df_3 %>%
      dplyr::filter(!!sym(color) %in% sams)
  }
  
  # 5' meta plots
  met_5 <- df_5 %>%
    create_meta(
      x         = "win_dist",
      y         = "counts",
      color     = color,
      alpha     = "type",
      plot_clrs = plot_clrs,
      plot_lvls = sams,
      se_clmn   = se_clmn,
      n_lab     = "n",
      ...
    ) +
    labs(y = plot_ttl) +
    add_breaks(seq(-10, 10), "TSS") +
    theme(
      legend.position    = "top",
      aspect.ratio       = 0.9,
      strip.text         = element_text(hjust = 0),
      panel.grid.major.x = element_line(linewidth = 0.5, color = "grey90")
    )
  
  # Remove legend if only one treatment
  if (n_distinct(df_5[[color]]) == 1) {
    met_5 <- met_5 +
      theme(legend.position = "none")
  }
  
  # 3' meta plots
  met_3 <- df_3 %>%
    create_meta(
      x         = "win_dist",
      y         = "counts",
      color     = "treat",
      alpha     = "type",
      plot_clrs = plot_clrs,
      plot_lvls = sams,
      se_clmn   = se_clmn,
      n_lab     = "n",
      ...
    ) +
    coord_cartesian(ylim = ylim_3) +
    add_breaks(seq(-3, 3, 1.5), "pAS") +
    theme(
      aspect.ratio       = 0.9,
      strip.text         = element_text(color = "white"),
      panel.grid.major.x = element_line(linewidth = 0.5, color = "grey90"),
      legend.position    = "none",
      axis.title.y       = element_blank()
    )
  
  # Facet wrap
  if (!is.null(grp)) {
    grp_f <- as.formula(str_c("~ ", grp))
    
    met_5 <- met_5 +
      facet_wrap(grp_f, ncol = 1)
    
    met_3 <- met_3 +
      facet_wrap(grp_f, ncol = 1)
  }
  
  # Return plot list
  if (return_list) {
    return(list(met_5, met_3))
  }
  
  # Create final figure
  res <- plot_grid(
    met_5, met_3,
    rel_widths = c(1, 1),
    nrow       = 1,
    align      = "vh",
    axis       = "trbl"
  )
  
  if (!is.null(file)) {
    ggsave(
      filename = basename(file),
      plot     = res,
      device   = dev,
      path     = dirname(file),
      width    = dims[1],
      height   = dims[2],
      units    = "in"
    )
  }
  
  res
}

# Helper to create 5' meta plot and pause stat panels
create_pausing_meta <- function(met_df, p_df, regions, sams, split_col = NULL, stats = c("fraction pause reads" = "pauses_NET"),
                                plot_clrs, file = NULL, dims = c(10, 6), v_line = c(0, 0.1, 0.3, 0.5, 1), bx_cols = plot_clrs,
                                bx_fill = "treat", med_pt = NULL, qnts = NULL) {
  
  # 5' metaplots
  p_met <- met_df %>%
    dplyr::filter(
      win_dist >= -1,
      treat %in% sams
    ) %>%
    create_meta(
      x          = "win_dist",
      y          = "counts",
      color      = "treat",
      alpha      = "type",
      plot_clrs  = plot_clrs,
      plot_lvls  = sams,
      v_line     = v_line,
      se_clmn    = "rep",
      vline_type = 2,
      n_lab      = "n"
    ) +
    labs(y = "mean signal (RPKM)") +
    add_breaks(seq(-10, 10), "TSS") +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme(
      strip.text      = element_text(hjust = 0),
      legend.position = "top",
      legend.text     = element_text(size = 8)
    )
  
  # Create boxplots showing pause stats
  bxs <- p_df %>%
    dplyr::filter(
      type   %in% stats,
      region %in% regions,
      treat  %in% sams
    ) %>%
    mutate(
      sample = fct_relevel(treat, names(plot_clrs)),
      type   = fct_relevel(type, stats)
    )
  
  bxs <- bxs %>%
    ggplot(aes(region, counts, fill = !!sym(bx_fill))) +
    geom_boxplot(alpha = 0.5, fatten = 1, outlier.size = 0.1, notch = TRUE) +
    
    geom_text(
      mapping       = aes(x = Inf, y = Inf, label = n),
      color         = "black",
      size          = txt_pt / .pt,
      alpha         = 1,
      check_overlap = T,
      hjust         = 1.2,
      vjust         = 1.4
    ) +
    
    facet_wrap(~ rep, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = bx_cols) +
    labs(y = names(stats)) +
    theme_info +
    theme(
      plot.title      = element_text(size = 10),
      legend.position = "none",
      axis.title.x    = element_blank(),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      axis.line       = element_blank(),
      panel.border    = element_rect(fill = NA, color = "black")
    )
  
  # Create final figure
  res <- plot_grid(
    p_met, bxs,
    rel_widths = c(0.75, 1),
    nrow       = 1,
    align      = "h",
    axis       = "tb"
  )
  
  if (!is.null(file)) {
    ggsave(
      filename = basename(file),
      plot     = res,
      device   = "png",
      path     = dirname(file),
      width    = dims[1],
      height   = dims[2],
      units    = "in"
    )
  }
  
  res
}

# Helper to create calculate pause stats and create pause stat panels
create_tss_pause_meta <- function(meta_df, pause_df, sams, meta_clrs, regions, bx_fill, bx_clrs = meta_clrs,
                                  key_vec = sam_lnms, metric = c("fraction pause reads" = "p-reads_NET")) {
  
  # Filter genes
  dat <- pause_df %>%
    dplyr::filter(
      sample %in% unname(sams),
      type == unname(metric)
    ) %>%
    
    get_overlapping_genes(c("sample", "region")) %>%  # remove genes that are not present in all samples and regions
    
    group_by(sample) %>%
    mutate(
      n = n_distinct(name),
      n = str_c("n = ", comma(n)),
    ) %>%
    ungroup() %>%
    
    mutate(
      region = recode(region, !!!regions),
      region = fct_relevel(region, regions),
      treat  = fct_relevel(treat, treats)
    )
  
  # Pause stats genes
  dat_genes <- dat %>%
    dplyr::select(name) %>%
    distinct()
  
  # Calculate mean mNET-seq signal
  p_mean_5 <- meta_df %>%
    semi_join(dat_genes, by = "name") %>%
    dplyr::filter(sample %in% names(sams))
  
  if (nrow(p_mean_5) == 0) {
    stop("No genes left after filtering, check input files")
  }
  
  p_mean_5 <- p_mean_5 %>%
    calc_mean_signal(
      rel  = FALSE,
      grps = c("sample", "group", "type")
    ) %>%
    format_sample_names(key_vec = key_vec)
  
  # Create final figure
  res <- create_pausing_meta(
    met_df        = p_mean_5,
    p_df          = dat,
    sams          = treats,
    regions       = unname(regions),
    stats         = metric,
    plot_clrs     = meta_clrs,
    bx_fill       = bx_fill,
    bx_cols       = bx_clrs
  )
  
  res
}

# Helper to create 5' zone meta plot and boxplot panels
create_zone_meta <- function(met_df, vln_df, sams, p_stat = "pauses_NET",  plot_cols,
                             file = NULL, dims = c(8, 5)) {
  
  # Create meta plots
  met <- met_df %>%
    dplyr::filter(
      type == p_stat,
      sample %in% sams
    ) %>%
    create_meta(
      plot_cols  = plot_cols,
      plot_lvls  = names(plot_cols),
      size       = 0.75,
      plot_title = NULL,
      h_line     = NULL,
      v_line     = NULL
    ) +
    
    geom_vline(
      aes(xintercept = zone, color = sample),
      linetype    = 2,
      show.legend = FALSE
    ) +
    
    scale_alpha_manual(values = 1) +
    
    labs(y = p_stat) +
    facet_wrap(~ group, scales = "free_y", ncol = 1) +
    theme(
      legend.position = "top",
      axis.title.x    = element_blank()
    )
  
  # Create violin plots
  vlns <- vln_df %>%
    dplyr::filter(sample %in% sams) %>%
    dplyr::select(group, sample, name, zone) %>%
    distinct() %>%
    group_by(group, sample) %>%
    mutate(n = n_distinct(name)) %>%
    ungroup() %>%
    mutate(sample = fct_relevel(sample, unique(names(plot_cols)))) %>%
    
    ggplot(aes(sample, zone, fill = sample)) +
    geom_boxplot(alpha = 0.5, fatten = 1, outlier.size = 0.2, notch = TRUE) +
    
    geom_text(
      mapping       = aes(x = Inf, y = Inf, label = str_c("n = ", n)),
      color         = "black",
      size          = txt_pt / .pt,
      alpha         = 1,
      check_overlap = T,
      hjust         = 1,
      vjust         = 1
    ) +
    
    facet_wrap(~ group, ncol = 1, scales = "free") +
    
    scale_fill_manual(values = plot_cols) +
    scale_y_log10() +
    labs(y = "pausing zone length (kb)") +
    
    theme_info +
    theme(
      legend.position = "none",
      strip.text      = element_text(color = "white"),
      axis.title.x    = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank()
    )
  
  # Create final figure
  res <- plot_grid(
    met, vlns,
    rel_widths = c(1, 0.4),
    align      = "h",
    axis       = "tb",
    nrow       = 1
  )
  
  if (!is.null(file)) {
    ggsave(
      filename = basename(file),
      plot     = res,
      device   = "png",
      path     = dirname(file),
      width    = dims[1],
      height   = dims[2],
      units    = "in"
    )
  }
  
  res
}


#' Filter data.frame and check for correctly formatted windows
#'
#' @param df_in data.frame containing counts for each gene window. Must at
#' least contain columns for gene name ('name'), counts ('counts'), and window
#' ID (specified by win_col).
#' @param win_num Expected number of windows for each gene. Set to NULL if a
#' variable number of windows is expected for each gene.
#' @param win_min Remove windows with an ID less than win_min.
#' @param win_max Remove windows with an ID greater than win_max.
#' @param win_len Expected length of each window. Set to NULL if different
#' lengths are expected for each window.
#' @param win_col Column containing window IDs.
#' @param count_lim Filter for genes where the sum of counts for all windows is
#' greater than count_lim.
#' @param filter_win_num Remove genes if they do not have the correct number of
#' windows specified by win_num.
#' @return Filtered data.frame
#' @export
.check_windows <- function(df_in, win_num = NULL, win_min = -Inf, win_max = Inf, win_len = NULL,
                           win_col = "win_id", count_lim = 0, filter_win_num = TRUE) {
  
  req_cols <- c("name", "counts", win_col)
  
  if (!all(req_cols %in% colnames(df_in))) {
    stop(str_c("Columns ", str_c(req_cols, collapse = ", "), " must be present in df_in."))
  }
  
  win_col <- sym(win_col)
  
  res <- df_in %>%
    group_by(name) %>%
    mutate(
      .nonsequential_wins = !all(min(!!win_col):max(!!win_col) %in% !!win_col),
      .incorrect_num_wins = FALSE,
      .incorrect_len_wins = FALSE
    )
  
  # Check genes for correct number of windows
  if (!is.null(win_num)) {
    res <- res %>%
      mutate(.incorrect_num_wins = n() != win_num)
  }
  
  if (filter_win_num && any(res$.incorrect_num_wins)) {
    n_removed <- res %>%
      dplyr::filter(.incorrect_num_wins) %>%
      pull(name) %>%
      n_distinct()
    
    warning(n_removed, " genes had the incorrect number of windows and were removed.")
    
    res <- res %>%
      dplyr::filter(!.incorrect_num_wins)
  }
  
  # Check genes for correct sized windows
  if (!is.null(win_len)) {
    if (!all(c("start", "end") %in% colnames(res))) {
      stop("data.frame must include 'start' and 'end' columns when win_len argument is provided.")
    }
    
    res <- res %>%
      mutate(.incorrect_len_wins = !all(end - start == win_len))
  }
  
  # Check for expected window format
  res %>%
    ungroup() %>%
    dplyr::select(starts_with(".")) %>%
    as.list() %>%
    iwalk(~ {
      if (any(.x)) {
        stop(str_c("Windows not formatted as expected: ", .y))
      }
    })
  
  # Filter windows
  res <- res %>%
    dplyr::filter(
      !!win_col >= win_min,
      !!win_col <= win_max
    ) %>%
    dplyr::filter(sum(counts) > count_lim) %>%
    ungroup() %>%
    dplyr::select(!starts_with("."))
  
  res
}

#' Calculate window distances
#'
#' @param df_in data.frame containing gene windows. Must contain columns for
#' gene name ('name') and window ID (specified by win_col).
#' @param win_len Length of windows. If NULL, df_in must include start and end
#' coordinates for each window.
#' @param ref_win ID for reference window. Distance values are calculated
#' relative to the ref_win. Upstream windows are reported as negative values.
#' @param win_col Column containing window IDs.
#' @param scale_dist Value to use for scaling distances. Set to 1000 to express
#' distances as kilobases.
#' @return data.frame with distances calculated for each gene window
#' @export
.calc_win_dist <- function(df_in, win_len = NULL, ref_win = 1, win_col = "win_id", scale_dist = 1000) {
  
  # Check data.frame columns
  req_cols <- c("name", win_col)
  
  coord_cols <- c("start", "end")
  
  if (is.null(win_len)) {
    req_cols <- c(coord_cols, req_cols)
  }
  
  if (!all(req_cols %in% colnames(df_in))) {
    stop(str_c("Columns ", str_c(req_cols, collapse = ", "), " must be present in df_in."))
  }
  
  # Check for window coordinates
  win_col <- sym(win_col)
  
  lens <- NULL
  
  coords_included <- all(coord_cols %in% colnames(df_in))
  
  if (coords_included) {
    lens <- unique(df_in$end - df_in$start)
  }
  
  # If win_len is provided
  # If window coordinates are also provided, check lengths
  if (!is.null(win_len) || length(lens) == 1) {
    if (!is.null(win_len) && coords_included) {
      stopifnot(
        length(lens) == 1,
        win_len == lens
      )
    }
    
    if (coords_included) {
      win_len <- lens
    }
    
    res <- df_in %>%
      mutate(
        win_dist = ((!!win_col - ref_win) * win_len) / scale_dist,
        win_len  = win_len / scale_dist
      ) %>%
      arrange(name, win_dist)
    
    return(res)
  }
  
  # If windows are different lengths calculate distance using cumsum
  res <- df_in %>%
    mutate(win_len = (end - start) / scale_dist) %>%
    arrange(name, !!win_col) %>%
    group_by(name) %>%
    mutate(
      id       = !!win_col,
      win_dist = cumsum(win_len),
      win_dist = lag(win_dist, 1, default = 0),
      ref_dist = ifelse(id == ref_win, win_dist, 0),
      ref_dist = max(ref_dist),
      win_dist = win_dist - ref_dist
    ) %>%
    ungroup() %>%
    dplyr::select(-ref_dist, -id)
  
  res
}

#' Merge bed files containing region windows so all files include the same genes
#' 
#' @param input_df data.frame or list of data.frames containing counts for each
#' gene window
#' @param groups If a data.frame is provided, columns to use for splitting
#' data.frame with dplyr::group_split().
#' @param win_num Expected number of windows for each gene. Set to NULL if a
#' variable number of windows is expected for each gene.
#' @param ref_win ID for reference window. Distance values are calculated
#' relative to the ref_win. Upstream windows are reported as negative values.
#' Set to NULL to skip this step.
#' @param win_min Remove windows with an ID less than win_min.
#' @param win_max Remove windows with an ID greater than win_max.
#' @param win_len Expected length of each window. Set to NULL if different
#' lengths are expected for each window.
#' @param win_col Column containing window IDs.
#' @param count_lim Filter for genes where the sum of counts for all windows is
#' greater than count_lim.
#' @param filter_win_num Remove genes if they do not have the correct number of
#' windows specified by win_num.
#' @param filter_unique Remove genes that are not shared between all samples.
#' @return Merged and filtered data.frame
#' @export
merge_wins <- function(input_df, groups = NULL, win_num = NULL, ref_win = 1, win_min = -Inf,
                       win_max = Inf, win_len = NULL, win_col = "win_id", count_lim = 0, 
                       filter_win_num = TRUE, filter_unique = TRUE, ...) {
  
  # Split data.frame into list based on groups
  if (is.data.frame(input_df)) {
    if (is.null(groups)) {
      stop("If data.frame is provided, must include groups for merging.")
    }
    
    input_df <- input_df %>%
      group_by(!!!syms(groups)) %>%
      group_split()
  }
  
  req_cols <- c("name", "counts", win_col)
  
  input_df %>%
    walk(~ {
      if (!all(req_cols %in% colnames(.x))) {
        req_cols <- str_c(req_cols, collapse = ", ")
        
        stop("Columns ", str_c(req_cols, collapse = ", "), " must be present in df_in.")
      }
    })
  
  # Filter and check windows
  res <- input_df %>%
    map(~ {
      .check_windows(
        df_in          = .x,
        win_num        = win_num,
        win_min        = win_min,
        win_max        = win_max,
        win_len        = win_len,
        win_col        = win_col,
        count_lim      = count_lim,
        filter_win_num = filter_win_num
      )
    })
  
  # Filter to only include genes present in all datasets
  if (filter_unique) {
    res <- res %>%
      merge_beds(...)
  }
  
  # Express window ID as distance based on a reference window
  if (!is.null(ref_win)) {
    res <- res %>%
      map(
        .calc_win_dist,
        ref_win    = ref_win,
        win_len    = win_len,
        win_col    = win_col,
        scale_dist = 1000
      )
  }
  
  res <- res %>%
    bind_rows() %>%
    dplyr::select(-any_of(c("chrom", "start", "end", "strand")))
  
  res
}

#' Merge bed files so all files include the same genes
#'
#' @param dfs List of data.frames to merge.
#' @param by Column containing containing gene names to use for merging.
#' @param bind Should data.frames be combined using bind_rows()? If set to
#' FALSE a list of data.frames will be returned.
#' @return List of data.frames only containing overlapping genes
#' @export
merge_beds <- function(dfs, by = "name", bind = FALSE) {
  
  # Create list of overlapping genes
  shared_genes <- dfs %>%
    map(distinct, !!sym(by)) %>%
    purrr::reduce(semi_join, by = by)
  
  if (length(shared_genes) == 0) {
    stop("There are no genes shared between the datasets.")
  }
  
  # Merge data.frames
  res <- dfs %>%
    map(semi_join, shared_genes, by = by)
  
  if (bind) {
    res <- res %>%
      bind_rows()
  }
  
  res
}

#' Run gprofiler
#' 
#' @param gene_list List of input genes.
#' @param genome Genome to use for identifying GO terms.
#' @param gmt_id GMT ID to use instead of genome name.
#' @param p_max Maximum adjusted p-value for GO term to be included in results.
#' @param GO_size Minimum term size for GO term.
#' @param intrsct_size Minimum number of genes that intersect GO term.
#' @param order_query Run an ordered query where more significant genes are at
#' the top of the list.
#' @param dbases GO databases to query.
#' @param file_path File path to write results.
#' @param overwrite Overwrite existing file.
#' @param ... Additional arguments to pass to gprofiler2.
#' @return tibble containing GO terms.
#' @export
run_gprofiler <- function(gene_list, genome = NULL, gmt_id = NULL, p_max = 0.05, GO_size = 10, intrsct_size = 10, 
                          order_query = FALSE, dbases = c("GO:BP", "GO:MF", "KEGG"), file_path = NULL,
                          overwrite = FALSE, ...) {
  
  # Check for empty gene list
  if (is_empty(gene_list)) {
    return(as_tibble(NULL))
  }
  
  # Check arguments
  if (is.null(genome) && is.null(gmt_id)) {
    stop("Must specifiy genome or gmt_id.")
  }
  
  if (!is.null(file_path) && file.exists(file_path) && !overwrite) {
    warning("Loading existing file: ", file_path, ". Set overwrite = TRUE to rerun.")
    
    return(read_tsv(file_path))
  }
  
  # Use gmt id
  if (!is.null(gmt_id)) {
    genome <- gmt_id
    dbases <- NULL
  }
  
  # Run gProfileR
  res <- gene_list %>%
    gost(
      organism       = genome,
      sources        = dbases,
      domain_scope   = "annotated",
      significant    = TRUE,
      user_threshold = p_max,
      ordered_query  = order_query,
      ...
    )
  
  # Format and sort output data.frame
  res <- as_tibble(res$result)
  
  if (!is_empty(res)) {
    
    if (!is.null(GO_size)) {
      res <- res %>%
        dplyr::filter(term_size > GO_size)
    }
    
    if (!is.null(intrsct_size)) {
      res <- res %>%
        dplyr::filter(intersection_size > intrsct_size)
    }
    
    res <- res %>%
      arrange(source, p_value)
  }
  
  # Write table
  if (!is.null(file_path) && nrow(res) > 0) {
    res %>%
      dplyr::select(-parents) %>%
      write_tsv(file_path)    
  }
  
  res
}

#' Create GO plot
#' 
#' @param GO_df data.frame of GO results.
#' @param plot_colors Plot colors.
#' @param n_terms Number of terms to label.
#' @return ggplot object
#' @export
create_bubbles <- function(go_df, plot_colors = NULL, n_terms = 15, txt_size = 6 / .pt, ...) {
  
  # Check for empty inputs
  if (is_empty(go_df) || nrow(go_df) == 0) {
    res <- ggplot() +
      geom_blank()
    
    return(res)
  }
  
  # Shorten GO terms and database names
  go_nms <- c(
    "GO:BP" = "Biological\nProcess",
    "GO:CC" = "Cellular\nComponent",
    "GO:MF" = "Molecular\nFunction",
    "KEGG"  = "KEGG"
  )
  
  go_dat <- go_df %>%
    mutate(
      term_id = str_remove(term_id, "(GO|KEGG):"),
      term_id = str_c(term_id, " ", term_name),
      term_id = str_to_lower(term_id),
      term_id = str_trunc(term_id, 30, "right"),
      source  = go_nms[source]
    )
  
  # Reorder database names
  plot_lvls <- unname(go_nms)
  plot_lvls <- plot_lvls[plot_lvls %in% go_dat$source]
  plot_lvls <- unique(c(plot_lvls, go_dat$source))
  
  go_dat <- go_dat %>%
    mutate(source = fct_relevel(source, plot_lvls))
  
  # Extract top terms for each database
  top_go <- go_dat %>%
    group_by(source) %>%
    arrange(p_value) %>%
    dplyr::slice(1:n_terms) %>%
    ungroup()
  
  # Create bubble plots
  go_cols <- c("#D55E00", "#0072B2", "#009E73", "#6A51A3")
  
  if (!is.null(plot_colors)) {
    go_cols <- plot_colors
  }
  
  res <- go_dat %>%
    ggplot(aes(1.25, -log10(p_value), size = intersection_size)) +
    geom_point(aes(color = source), alpha = 0.5, show.legend = TRUE, ...) +
    
    geom_text_repel(
      aes(2, -log10(p_value), label = term_id),
      data         = top_go,
      size         = txt_size,
      direction    = "y",
      hjust        = 0,
      segment.size = NA
    ) +
    
    guides(color = "none") +
    xlim(1, 8) +
    labs(y = "-log10(p-value)") +
    
    scale_color_manual(values = go_cols) +
    theme_info +
    theme(
      legend.position = "bottom",
      axis.title.x    = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.line.x     = element_blank()
    ) +
    facet_wrap(~ source, scales = "free", nrow = 1)
  
  res
}

# Function to create sequence logos
create_logo <- function(input_freq, plot_colors, x_lim = NULL, y_lim = NULL,
                        x_breaks = NULL, x_labs = NULL, v_line = NULL) {
  
  # Set color scheme
  plot_cols <- make_col_scheme(
    chars = c("A", "T", "G", "C"), 
    cols  = plot_colors
  )
  
  # Simplify table 
  freq_df <- input_freq %>% 
    mutate(
      A    = round(A, 3), 
      `T`  = round(`T`, 3), 
      G    = round(G, 3), 
      C    = round(C, 3),
      tot  = (A + `T` + G + C),
      diff = 1 - tot,
      C    = ifelse(tot != 1, C + diff, C)
    ) %>%
    dplyr::select(-tot, -diff) %>%
    gather(nuc, value, -position)
  
  # Set x-limits 
  if (!is.null(x_lim)) {
    x_min <- x_lim[1]
    x_max <- x_lim[2]
    
    freq_df <- freq_df %>% 
      dplyr::filter(position >= x_min & position <= x_max)
  }
  
  freq_df <- freq_df %>% 
    spread(position, value) 
  
  # Add rownames
  rownames(freq_df) <- freq_df$nuc
  
  freq_df <- freq_df %>% 
    dplyr::select(-nuc)
  
  # Create sequence logo
  res <- ggseqlogo(as.matrix(freq_df), col_scheme = plot_cols) +
    scale_y_continuous(breaks = seq(0, 2, 0.5))
  
  # Change x-axis labels
  if (!is.null(x_breaks)) {
    x_breaks <- seq(1 : (x_lim[2] - x_lim[1] + 1))
    
    res <- res +
      scale_x_continuous(
        breaks = x_breaks, 
        labels = x_labs
      )
  }
  
  # Add vertical dotted line
  if (!is.null(v_line)) {
    res <- res + 
      geom_vline(
        xintercept = v_line, 
        linetype   = 2
      )
  }
  
  # Set y-limits
  if (!is.null(y_lim)) {
    res <- res + 
      coord_cartesian(ylim = y_lim)
  }
  
  res
}

# Calculate pause stats
calc_pause_stats <- function(df_in, n_col = "NET", r_col = "p-reads", p_col = "pauses",
                             len_col = "len", repl_na = NULL, log_trans = NULL, len_div = 1000) {
  
  stat_cols  <- c(n_col, r_col, p_col)
  other_cols <- colnames(df_in)
  other_cols <- other_cols[!other_cols %in% stat_cols]
  
  # Add pausing stats
  # NAs will be generated for regions with no NET-seq reads
  res <- df_in %>%
    mutate(
      across(
        all_of(c(n_col, r_col, p_col)),     # normalized by length
        ~ .x / (!!sym(len_col) / len_div),
        .names = "{.col}_kb"
      ),
      across(
        all_of(c(p_col, r_col)),            # normalized by NET reads
        ~ .x / !!sym(n_col),
        .names = "{.col}_{n_col}"
      ),
      across(
        all_of(r_col),                      # normalized by pauses
        ~ .x / !!sym(p_col),
        .names = "{.col}_{p_col}"
      )
    )
  
  # Replace NAs
  # could observe NAs in pauses_NET and p-reads_NET if 0 NET-seq reads
  # could observe NAs in p-reads_pauses if 0 pauses
  if (!is.null(repl_na)) {
    res <- res %>%
      mutate(across(
        all_of(repl_na),
        replace_na,
        0
      ))
  }
  
  # Log transform values
  if (!is.null(log_trans)) {
    res <- res %>%
      mutate(across(
        all_of(log_trans),
        log10,
        .names = "log10_{.col}"
      ))
  }
  
  # Format output data.frame
  res <- res %>%
    pivot_longer(
      cols      = -all_of(other_cols),
      names_to  = "type",
      values_to = "counts"
    )
  
  res
}

# Function to add subscripts to plot labels
add_subscript <- function(input_str, pattern = "DN|WT", bold = T) {
  
  # Extract characters and position of pattern
  pos     <- str_locate(input_str, pattern)
  str_len <- str_length(input_str)
  chars   <- str_split(input_str, "") %>% 
    unlist()
  
  # Convert pattern to subscript
  if (str_detect(input_str, pattern)) {
    start_chars <- ""
    pattern     <- chars[pos[1] : pos[2]] %>% 
      str_c(collapse = "")
    
    if (pos[1] > 1) {
      start_chars <- chars[1 : (pos[1] - 1)] %>%
        str_c(collapse = "")
    }
    
    end_chars <- ""
    
    if (pos[2] < str_len) {
      end_chars <- chars[(pos[2] + 1) : str_len] %>%
        str_c(collapse = "")
    }
    
    res <- bquote(.(start_chars) * ""[.(pattern)] * .(end_chars))
    
    if (bold) {
      res <- bquote(bold(.(start_chars) * ""[.(pattern)] * .(end_chars)))
    }
    
  } else {
    start_chars <- chars[1 : str_len] %>%  # This was only way to get facet_wrap labeller working
      str_c(collapse = "")
    
    res <- start_chars
    
    if (bold) {
      res <- bquote(bold(.(start_chars)))
    }
  }
  
  res
}

# Cluster genes
cluster_genes <- function(df_in, gene_col = "name", k = 8, n_top = NULL,
                          seed = 42, ...) {
  set.seed(seed)
  
  clust_genes <- df_in %>%
    pull(gene_col)
  
  clusts <- df_in %>%
    as.data.frame() %>%
    dplyr::select(-!!sym(gene_col)) %>%
    kmeans(centers = k, ...)
  
  clusts <- clusts$cluster %>%
    tibble(
      name  = clust_genes,
      clust = .
    )
  
  # Select top clusters and merge others
  if (!is.null(n_top)) {
    top_clusts <- clusts %>%
      group_by(clust) %>%
      summarize(n = n_distinct(name), .groups = "drop") %>%
      slice_max(n, n = n_top) %>%
      pull(clust) %>%
      set_names(1:n_top, .)
    
    clusts <- clusts %>%
      mutate(
        top = clust %in% names(top_clusts),
        clust = recode(clust, !!!top_clusts),
        clust = ifelse(top, clust, 3)
      ) %>%
      dplyr::select(-top)
  }
  
  clusts
}

# Set y-axis scales so signal ratio for strands is equal for metaplots
set_strand_ratio <- function(df_in, df_ref, mtplyr = 3, count_col = "counts") {
  
  # Calculate sense/anti-sense ratio
  df_ref <- df_ref %>%
    mutate(strand = if_else(!!sym(count_col) >= 0, "pos", "neg"))
  
  st_max <- df_ref %>%
    group_by(strand) %>%
    summarize(max = max(abs(!!sym(count_col))), .groups = "drop")
  
  st_max <- set_names(st_max$max, st_max$strand)
  st_rat <- max(st_max) / min(st_max)
  
  # Determine scale range
  sig_rng <- max(pull(df_in, count_col)) * mtplyr
  sig_rng <- c(sig_rng / -st_rat, sig_rng)
  
  # Determine which strand has higher signal
  st_top <- st_max[st_max == max(st_max)] %>%
    names()
  
  if (st_top == "neg") {
    sig_rng <- rev(sig_rng * -1)
  }
  
  sig_rng
}


#' Create browser plots
#' 
#' @param bgs List of data.frame lists containing bedGraphs for each sample
#' for each strand
#' @param genes_df data.frame containing gene coordinates to use for gene of
#' interest
#' @param genes_db EnsDb containing exon info for each gene
#' @param track_clrs Colors to use for each sample
#' @param strand Strand to plot, if NULL signal for the strand matching the
#' provided gene will be plottted.
#' @param pauses List of data.frames containing pause coordinates for each
#' sample
#' @param zones List of data.frames containing pausing zone coordinates for
#' each sample
#' @param ln_clr Line color to use for marking pausing zone
#' @param track_scale Fraction of region to show for track scale
#' @param flank Vector containing the fraction of gene length to use for
#' expanding plot coordinates
#' @param track_ttls Vector containing titles to label tracks
#' @param equal_axis Set y-axis limits equal for all tracks
#' @param y_dec When equal_axis is TRUE, number of decimals to round axis
#' labels
#' @param pause_track_height Relative height of the pause annotation track
#' @param gene_track_height Relative height of the gene annotation track
#' @param scale_track_height Relative height of the gene scale track
#' @param line_width Width of line for histogram
create_browser <- function(bgs, gene, genes_df, genes_db, track_clrs,
                           strand = NULL, pauses = NULL, zones = NULL,
                           ln_clr = "#E62D30", track_scale = 0.1,
                           flank = c(0.05, 0.05), track_ttls = NULL,
                           equal_axis = FALSE, y_dec = 1,
                           scale_track_height = 0.1, pause_track_height = 0.15,
                           gene_track_height = 0.15, line_width = 2, ...) {
  
  # Retrieve gene info
  # fix chromosome names
  goi <- genes_df %>%
    dplyr::filter(grepl(str_c("(?<=\\|)", gene, "$"), name, perl = TRUE)) %>%
    mutate(gene_len = end - start)
  
  if (nrow(goi) != 1) {
    stop("gene must return exactly one row from genes_df")
  }
  
  chr <- str_remove(goi$chrom, "^chr")
  
  # Retrieve exons to use for gene model
  # use exon set from transcript with closest coordinates
  tx_id <- genes_db %>%
    transcripts(filter = GeneNameFilter(gene)) %>%
    as_tibble() %>%
    mutate(
      start_dif = abs(start - goi$start),
      coord_dif = start_dif + abs(end - goi$end)
    ) %>%
    dplyr::filter(coord_dif == min(coord_dif)) %>%
    dplyr::filter(start_dif == min(start_dif)) %>%
    dplyr::filter(tx_id == dplyr::first(tx_id)) %>%
    pull(tx_id)
  
  goi_track <- genes_db %>%
    getGeneRegionTrackForGviz(filter = TxNameFilter(tx_id))
  
  # Gene models
  mod_clr <- "#09004C"
  
  grtrack <- GeneRegionTrack(
    range      = goi_track,
    genome     = "hg38",
    chromosome = chr,
    gene       = gene,
    col        = mod_clr,
    fill       = mod_clr,
    fontcolor  = mod_clr,
    col.line   = mod_clr
  )
  
  # Axis scale
  axisTrack <- GenomeAxisTrack(chromosome = goi$chrom, col = "black")
  
  # Bigwig data
  # fix chromosome names
  track_dat <- bgs[[goi$strand]]
  
  track_dat <- track_dat %>%
    map(~ {
      .x %>%
        dplyr::filter(chrom == goi$chrom) %>%
        dplyr::select(-chrom)
    })
  
  flank[1] <- round(goi$gene_len * flank[1])
  flank[2] <- round(goi$gene_len * flank[2])
  
  track_heights <- rep(scale_track_height, rep(1, length(track_dat)), gene_track_height)
  
  # Get max signal for track region
  region <- goi$start - flank[1]
  region <- c(region, goi$end + flank[2])
  
  y_ticks <- track_dat %>%
    map_dbl(~ {
      val <- .x %>%
        dplyr::filter(
          start >= region[1],
          end   <= region[2]
        ) %>%
        pull(value)
      
      ifelse(any(val >= 0), max(val), min(val))
    })
  
  if (equal_axis) {
    if (any(y_ticks >= 0)) {
      y_ticks <- c(0, max(y_ticks))
    } else {
      y_ticks <- c(min(y_ticks), 0)
    }
    
    y_lim <- vector("list", length(track_dat))
    
    y_lim[seq_along(y_lim)] <- list(round(y_ticks, y_dec))
    
    y_ticks <- y_lim
    
  } else {
    y_lim <- y_ticks <- y_ticks %>%
      map(~ {
        val <- round(.x, y_dec)
        
        if (val >= 0) return(c(0, val))
        
        c(val, 0)
      })
  }
  
  # Create signal tracks
  track_args <- list(
    dat   = track_dat,
    clrs  = track_clrs[names(track_dat)],
    nms   = track_ttls %||% names(track_dat),
    lim   = y_lim,
    ticks = y_ticks
  )
  
  dtracks <- track_args %>%
    pmap(~ {
      args <- list(...)
      
      DataTrack(
        range          = args$dat,
        genome         = "hg38",
        chromosome     = chr,
        strand         = goi$strand,
        type           = "histogram",
        windowSize     = 1,
        lwd            = line_width,
        fill.histogram = args$clrs,
        col.histogram  = args$clrs,
        name           = args$nms,
        ylim           = args$lim,
        yTicksAt       = args$ticks
      )
    })
  
  # Highlight pause zones
  if (!is.null(zones)) {
    dtracks <- dtracks %>%
      imap(~ {
        starts <- zones[[.y]] %>%
          dplyr::filter(grepl(str_c("(?<=\\|)", gene, "$"), name, perl = TRUE)) %>%
          pull(zone)
        
        if (is_empty(starts)) {
          return(.x)
        }
        
        HighlightTrack(
          trackList  = .x,
          genome     = "hg38",
          chromosome = chr,
          start      = starts,
          width      = 0,
          lwd        = 2,
          lty        = 2,
          col        = ln_clr,
          fill       = ln_clr
        )
      })
  }
  
  # Highlight pause sites
  if (!is.null(pauses)) {
    dtracks <- dtracks %>%
      imap(~ {
        starts <- pauses[[.y]] %>%
          dplyr::filter(grepl(str_c("(?<=\\|)", gene, "$"), name, perl = TRUE)) %>%
          pull(start)
        
        ann_track <- AnnotationTrack(
          genome     = "hg38",
          chromosome = chr,
          strand     = goi$strand,
          start      = starts,
          width      = 0,
          stacking   = "dense",
          col        = NA,
          shape      = "box"
        )
        
        list(.x, ann_track)
      })
    
    dtracks <- flatten(dtracks)
    
    track_heights <- c(
      scale_track_height,
      rep(c(1, pause_track_height), length(track_dat)),
      gene_track_height
    )
  }
  
  # Create figure
  tracklst <- list(axisTrack) %>%
    append(dtracks) %>%
    append(list(grtrack))
  
  plotTracks(
    trackList        = tracklst,
    from             = goi$start,
    to               = goi$end,
    extend.left      = flank[1],
    extend.right     = flank[2],
    background.title = NA,
    fontcolor.title  = "black",
    fontface.title   = "plain",
    col.axis         = "black",
    face.axis        = "plain",
    scale            = track_scale,
    sizes            = track_heights,
    ...
  )
}

.get_bg_clrs <- function(samples, plt_grps, sample_df) {
  
  res <- sample_df %>%
    dplyr::filter(
      sample %in% samples,
      plot_grp %in% plt_grps
    )
  
  res <- set_names(
    res$clr, res$sample
  )
  
  res
}

.get_bgs <- function(samples, grps, plt_grps, sample_df,
                     subsample_region = "_gene", results_dir = res_dir) {
  
  args <- list(
    samples  = samples,
    grps     = grps,
    plt_grps = plt_grps
  )
  
  bg_paths <- args %>%
    pmap_chr(~ {
      args <- list(...)
      
      fl <- sample_df %>%
        dplyr::filter(
          sample       == args$samples,
          sampling_grp == args$grps,
          plot_grp     == args$plt_grps
        ) %>%
        mutate(file = str_c(file, "-", sampling_grp, subsample_region)) %>%
        pull(file)
      
      file.path(results_dir, fl, fl)
    })
  
  names(bg_paths) <- args$samples
  
  strnds <- list("+" = "pos", "-" = "neg")
  
  res <- strnds %>%
    map(~ {
      strd <- .x
      
      bg_paths %>%
        map(~ read_bedgraph(str_c(.x, "_", strd, "_N.bedgraph.gz"))) %>%
        map(mutate, start = start + 1)
    })
  
  res
}

.get_pauses <- function(samples, grps, plt_grps, sample_df,
                        subsample_region = "_gene", results_dir = res_dir) {
  
  args <- list(
    samples  = samples,
    grps     = grps,
    plt_grps = plt_grps
  )
  
  paths <- args %>%
    pmap_chr(~ {
      args <- list(...)
      
      sample_df %>%
        dplyr::filter(
          sample       == args$samples,
          sampling_grp == args$grps,
          plot_grp     == args$plt_grps
        ) %>%
        mutate(file = str_c(file, "-", sampling_grp, subsample_region)) %>%
        pull(file)
    })
  
  names(paths) <- samples
  
  # must convert to 1-based coordinates
  res <- paths %>%
    map_chr(~ here(results_dir, .x, "pauses", str_c(.x, "_200_strong_pauses.bed.gz")))
  
  res <- res %>%
    map(read_bed, n_fields = 6) %>%
    map(mutate, start = start + 1)
  
  res
}

.get_zones <- function(plt_grps, genes_df, results_dir = plot_params$obj_dir,
                       prefix = str_c("_", plot_params$pause_win, plot_params$pause_strength)) {
  
  # Load pausing zones
  res <- unique(plt_grps) %>%
    map_dfr(~ {
      fl <- str_c(.x, prefix, "pause_zones.tsv.gz")
      fl <- here(results_dir, fl)
      
      fl %>%
        vroom() %>%
        distinct(
          name,  group, sample,
          treat, rep,   zone,
          zone_class
        )
    })
  
  res <- genes_df %>%
    dplyr::select(chrom, start, end, name, strand) %>%
    inner_join(res, by = "name")
  
  # Set pausing zone coordinates
  # must convert to 1-based coordinates
  res <- res %>%
    dplyr::filter(zone_class == "zone") %>%
    mutate(
      zone_len = zone,
      zone     = zone * 1000,
      zone     = ifelse(strand == "+", start + 1 + zone, end - zone)
    )
  
  res <- res %>%
    split(.$sample)
  
  res
}


# Helper to load pause bed files
load_pause_beds <- function(prfxs, sfxs, paths, genes = NULL, regions = NULL) {
  # paths <- str_c(paths, "/", prfxs, sfxs)
  # paths <- set_names(paths, prfxs)
  
  res <- paths %>%
    map(read_tsv, col_names = bed_cols)
  
  res <- res %>%
    imap_dfr(~ {
      .x <- .x %>%
        separate(
          name,
          into   = c("pause_coords", "gene_name"),
          sep    = "\\*",
          remove = FALSE
        ) %>%
        mutate(
          sample = .y,
          name   = paste0(name, "*", score)
        ) %>%
        dplyr::select(-pause_coords)
      
      if (!is.null(genes)) {
        .x <- .x %>%
          semi_join(genes, by = c("gene_name" = "name"))
      }
      
      .x
      
      # THIS PRODUCES SPORADIC ERROR WHEN KNITTING: "embedded nul in string"
      # .x %>%
      #   mutate(
      #     sample    = .y,
      #     gene_name = str_extract(name, "(?<=\\*).+$"),
      #     name      = paste0(name, "*", score)
      #   ) %>%
      #   semi_join(genes, by = c("gene_name" = "name"))
    })
  
  n_pauses <- nrow(res)
  
  # Intersect pauses
  if (!is.null(regions)) {
    res <- regions %>%
      imap_dfr(~ {
        bed_intersect(res, .x, suffix = c("", ".y")) %>%
          dplyr::filter(
            .overlap == 1,      # by default bookended intervals are included
            strand == strand.y  # some pAS regions from opposite strands will overlap between genes
          ) %>%
          group_by(sample) %>%
          mutate(
            region = .y,
            n_lab  = str_c(comma(n_distinct(name)), " pauses\n"),
            n_lab  = str_c(n_lab, comma(n_distinct(gene_name)), " genes"),
          ) %>%
          dplyr::select(-ends_with(".y"), -.overlap)
      }) %>%
      ungroup()
    
    if (n_pauses != nrow(res)) {
      stop(
        "The number of pauses after intersection does not match the number of ",
        "input pauses, data.frame is malformed. This could be due to overlap ",
        "between the provided regions or overlap between genes."
      )
    }
  }
  
  res
}

# Set sample colors
# basic test for valid rgb color and type and
# switches hex <--> rgb, able to apply tint colors
RgbToHex <- function(x, 
                     convert = "hex", 
                     tint = FALSE) {
  if(str_count(x, ",") == 2){
    # hex <- rgb
    myhex <- try(rgb(str_split_fixed(x,",",n=3),
                     maxColorValue = 255), silent = TRUE)
    if("try-error" %in% class(myhex)){
      myhex <- "#000000"
      mygrb <- "0,0,0"
    } else {
      myrgb <- x
    }
  } else {
    myrgb <- try(col2rgb(x), silent = TRUE)
    if("try-error" %in% class(myrgb)){
      myhex <- "#000000"
      myrgb <- "0,0,0"
    } else {
      myhex <- x
      myrgb <- paste(myrgb, collapse = ",")
    }
  }
  if (is.numeric(tint)) {
    myhex <- lighten(myhex,amount = 0.2 * (tint-1))
    myrgb <- col2rgb(myhex)
  } 
  if(convert == "hex"){
    return(myhex)
  } else {
    return(myrgb) 
  }
}

# Subsample sample directories
get_file_dir <- function(sam_dir, file_dir) {
  sam_dir %>%
    set_names() %>%
    map_chr(~ here(res_dir, file_dir))
}

# Gene lists and filter based on length
load_genes <- function(drctry, file, len = 1, rm_cols = TRUE) {
  res <- here(drctry, file) %>%
    read_tsv(col_names = bed_cols) %>%
    mutate(gene_len = (end - start) / 1000) %>%
    dplyr::filter(gene_len > len)
  
  if (rm_cols) {
    res <- res %>%
      dplyr::select(name, gene_len)
  }
  
  res
}

# Add extra values so signal appears as bar graph
convert_to_bars <- function(df_in, xclmn = "win_id", yclmn = "counts") {
  res <- df_in %>%
    rowwise() %>%
    mutate(
      !!sym(xclmn) := list(c(!!sym(xclmn) - 0.5, !!sym(xclmn) + 0.5)),
      !!sym(yclmn) := list(rep(!!sym(yclmn), 2))
    ) %>%
    ungroup() %>%
    unnest(all_of(c(xclmn, yclmn)))
  
  res
}

# Calculate overlap between pause bed files
calc_pause_overlap <- function(dfs) {
  
  if (length(dfs) < 2) {
    return(NULL)
  }
  
  # Identify shared pauses
  y_i <- 2:length(dfs)
  
  shared <- bed_intersect(
    x      = dfs[[1]],
    y      = dfs[y_i],
    suffix = c("", ".y")
  ) %>%
    dplyr::filter(
      strand == strand.y,
      .overlap > 0
    )
  
  # If more than two samples get pauses shared between all
  if (length(dfs) > 2) {
    n_y <- seq_along(y_i)
    
    shared <- shared %>%
      group_by(name) %>%
      dplyr::filter(all(n_y %in% .source)) %>%
      ungroup()
  }
  
  shared <- shared %>%
    dplyr::select(-ends_with(".y"), -starts_with(".")) %>%
    distinct()
  
  n_shared <- nrow(shared)
  
  res <- tibble(
    key      = str_c(names(dfs), collapse = " "),
    class    = "shared",
    n_pauses = n_shared,
    n_genes  = n_distinct(shared$gene_name),
    n        = n_pauses
  )
  
  # Identify unique pauses
  uniq <- tibble(
    key      = names(dfs),
    class    = "unique",
    n_pauses = map_int(dfs, nrow),
    n_genes  = map_int(dfs, ~ n_distinct(.x$gene_name)),
    n        = n_pauses - n_shared
  )
  
  # Create results table
  res <- bind_rows(res, uniq) %>%
    mutate(frac = n / sum(n))
  
  res
}

# Calculate mean signal for pause metaplots
calc_mean_pause_signal <- function(df_in, grps, key_vec = sam_lnms, rel = FALSE) {
  res <- df_in %>%
    group_by(!!!syms(grps)) %>%
    mutate(
      orig_n = n_lab,
      n_lab  = str_c(comma(n_distinct(name)), " pauses\n"),
      n_lab  = str_c(n_lab, comma(n_distinct(gene_name)), " genes"),
    ) %>%
    ungroup()
  
  if (!identical(res$n_lab, res$orig_n)) {
    warning("The number of genes and/or pauses per group is different from the original calculation")
  }
  
  res <- res %>%  
    calc_mean_signal(
      win_col = "win_id",
      rel     = rel,
      grps    = c(grps, "n_lab")
    ) %>%
    format_sample_names(key_vec = key_vec) %>%
    convert_to_bars() %>%
    mutate(win_dist = win_id - 50)
  
  res
}

# Plot mean signal
create_pause_meta <- function(df_in, x, y, color = NULL, plot_clrs = NULL, size = 1, v_line = 0,
                              vline_clr = "black", vline_type = 2, n_lab = NULL, n_lab_pos = c(Inf, Inf),
                              x_lab = waiver(), y_lab = waiver(), ...) {
  res <- df_in
  
  if (!is.null(names(plot_clrs)) && !is.null(color)) {
    res <- res %>%
      mutate(!!sym(color) := fct_relevel(!!sym(color), names(plot_clrs)))
  }
  
  res <- res %>%
    ggplot(aes(!!sym(x), !!sym(y), color = !!sym(color), fill = !!sym(color))) +
    geom_ribbon(
      aes(ymin = 0, ymax = !!sym(y)),
      alpha        = 0.2,
      linewidth    = 0.5,
      outline.type = "upper",
      key_glyph    = draw_key_point,
      ...
    ) +
    guides(fill = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    geom_vline(xintercept = v_line, linetype = vline_type, color = vline_clr) +
    scale_color_manual(values = plot_clrs) +
    scale_fill_manual(values = plot_clrs) +
    theme_info +
    theme(
      legend.position = "none",
      axis.line       = element_blank(),
      panel.border    = element_rect(fill = NA, color = "black")
    ) +
    labs(x = x_lab, y = y_lab)
  
  if (!is.null(n_lab)) {
    res <- res +
      geom_text(
        mapping       = aes(x = n_lab_pos[1], y = n_lab_pos[2], label = !!sym(n_lab)),
        color         = "black",
        size          = txt_pt / 1.2 / .pt,
        alpha         = 1,
        check_overlap = TRUE,
        hjust         = 1.1,
        vjust         = 1.15
      )
  }
  
  res
}

# Create pause logo
create_pause_logo <- function(df_in, left = 9, right = 10, genome, vline = 10.5) {
  
  cen <- vline - 0.5
  
  n_lab <- unique(df_in$n_lab)
  
  res <- df_in %>%
    bed_slop(left = left, right = right, genome = genome, strand = TRUE) %>%
    dplyr::rename(chr = chrom) %>%
    dplyr::select(chr, start, end, name, score, strand)
  
  res <- res %>% 
    as.data.frame() %>%
    get.fasta(
      fasta   = fa,
      strand  = TRUE,  # must include all 6 columns for strand option to work
      verbose = FALSE
    ) %>%
    .$sequence
  
  res <- res %>%
    ggseqlogo(seq_type = "dna") +
    
    geom_vline(xintercept = vline, linetype = 2) +
    
    scale_x_continuous(
      labels = (function(x) ifelse(x <= cen, x - (cen + 1), str_c("+", x - cen))),
      breaks = c(1, 6, 15, 20)
    ) +
    
    scale_y_continuous(breaks = breaks_limits) +
    theme_info +
    theme(
      axis.line = element_blank(),
      panel.border = element_rect(fill = NA, color = "black")
    )
  
  res
}

# Create figure showing pause signal and pause logo
create_pause_logo_fig <- function(met_df, pause_df, ttl, y_lab, y_lim = NULL, genome) {
  
  met <- met_df %>%
    create_pause_meta(
      x         = "win_dist",
      y         = "counts",
      color     = "treat",
      v_line    = 0.5,
      plot_clrs = clrs,
      n_lab     = "n_lab",
      x_lab     = "distance from pause (bp)",
      y_lab     = y_lab
    ) +
    coord_cartesian(ylim = y_lim) +
    scale_y_continuous(breaks = breaks_limits) +
    ggtitle(ttl) +
    theme(
      plot.title   = element_text(size = ttl_pt),
      plot.margin  = margin(7, 7, 0, 7, "pt"),
      axis.ticks.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.title.x = element_blank()
    )
  
  lo <- pause_df %>%
    create_pause_logo(genome = genome) +
    theme(plot.margin = margin(0, 7, 7, 7, "pt"))
  
  res <- plot_grid(
    met, lo,
    nrow        = 2,
    align       = "v",
    rel_heights = c(1, 0.5)
  )
  
  res
}

get_box_data <- function(df_in, p_lim, tss_lim, body_lim, only_shared = FALSE) {
  
  # Filter based on number of pauses in each region
  mets <- unique(df_in$type)
  
  res <- df_in %>%
    pivot_wider(names_from = "type", values_from = "counts") %>%
    dplyr::filter(
      pauses >= p_lim,
      NET    >= ifelse(region == "body", 1, tss_lim),
      NET_kb >= ifelse(region == "body", body_lim, 0)
    ) %>%
    pivot_longer(all_of(mets), names_to = "type", values_to = "counts")
  
  # Remove genes that do not have data for all regions
  # do not remove rows with NaNs since these occur for p-reads_pauses when
  # there are zero pauses
  if (only_shared) {
    regs <- unique(res$region)
    
    res <- res %>%
      pivot_wider(names_from = "region", values_from = "counts") %>%
      
      dplyr::filter(if_all(all_of(regs), ~ is.nan(.x) | !is.na(.x))) %>%
      
      pivot_longer(all_of(regs), names_to = "region", values_to = "counts")
  }
  
  # Remove regions that do not have data for all treatment groups
  # do not remove rows with NaNs since these occur for p-reads_pauses when
  # there are zero pauses
  treats <- unique(res$treat)
  
  res <- res %>%
    dplyr::select(-sample) %>%
    pivot_wider(names_from = "treat", values_from = "counts") %>%
    
    dplyr::filter(if_all(all_of(treats), ~ is.nan(.x) | !is.na(.x))) %>%
    
    pivot_longer(all_of(treats), names_to = "treat", values_to = "counts")
  
  res
}

create_stat_boxes <- function(df_in, metric, by = c("rep", "type", "region"),
                              plt_clrs = clrs, y_trans = "identity",
                              reg_labels = regs, p_treats = fc_treats, lvls = treats) {
  
  reg_labels <- reg_labels[names(reg_labels) %in% df_in$region]
  
  treats <- unique(df_in$treat)
  n_trts <- length(treats)
  
  # Check that same genes are present for all samples  
  dat <- df_in %>%
    dplyr::filter(type == metric) %>%
    pivot_wider(names_from = "treat", values_from = "counts")
  
  if (any(!complete.cases(dat))) {
    stop("NAs present in data.")
  }
  
  # Format data for boxplots
  # add n labels
  dat <- dat %>%
    group_by(!!!syms(by)) %>%
    mutate(
      n      = n_distinct(name),
      n      = str_c("n = ", comma(n)),
      bx_lab = n
    )
  
  n_n <- n_distinct(dat$n)
  
  # Calculate p-value for treats
  if (n_trts > 1) {
    n_tests <- n_groups(dat)
    
    dat <- dat %>%
      mutate(
        p = list(broom::tidy(
          wilcox.test(!!!syms(p_treats))
        )),
        p = map_dbl(p, ~ .x$p.value)
      )
    
    # Correct for multiple testing
    dat <- dat %>%
      rowwise() %>%
      mutate(
        p = p.adjust(p, method = "bonferroni", n = n_tests),
        p = scientific(p, digits = 2),
        bx_lab = p
      ) %>%
      ungroup()
    
    if (n_n > 1) {
      dat <- dat %>%
        mutate(bx_lab = str_c(n, "\n", p))
    }
  }
  
  dat <- dat %>%
    ungroup() %>%
    pivot_longer(
      all_of(treats),
      names_to  = "treat",
      values_to = "counts"
    ) %>%
    mutate(
      treat  = fct_relevel(treat, lvls),
      region = reg_labels[region],
      region = fct_relevel(region, unname(reg_labels))
    )
  
  # Create boxplots
  res <- dat %>%
    ggplot(aes(region, counts, fill = treat)) +
    
    geom_boxplot(
      alpha        = 0.5,
      fatten       = 1,
      outlier.size = 0.1,
      notch        = TRUE,
      key_glyph    = draw_key_point
    ) +
    
    scale_fill_manual(values = plt_clrs) +
    
    guides(fill = guide_legend(override.aes = list(shape = 22, size = 4, stroke = 1, color = "black"))) +
    labs(y = names(metric)) +
    box_theme
  
  # Add n labels
  y_exp  <- 0.05
  v_jst  <- 0
  lab_pt <- txt_pt / .pt
  
  if (n_trts > 1 || n_n > 1) {
    y_exp <- y_exp + 0.15
    
    v_jst <- v_jst + if_else(n_n > 1, 1.2, 1.5)
    
    res <- res +
      geom_text(
        mapping       = aes(x = region, y = Inf, label = bx_lab),
        color         = "black",
        size          = lab_pt,
        alpha         = 1,
        check_overlap = TRUE,
        vjust         = v_jst
      )
    
    v_jst <- v_jst + 1.5
  }
  
  if (identical(n_n, 1L)) {
    # n_lab <- unique(dat$n)
    # res <- res +
    #   labs(subtitle = n_lab) +
    #   theme(plot.subtitle = element_text(size = lab_pt * .pt, hjust = 0.98, vjust = 0))
    
    if (identical(n_trts, 1L)) {
      lab_pt <- lab_pt * 1.5
    }
    
    y_exp <- y_exp + 0.1
    v_jst <- v_jst + 1.5
    
    res <- res +
      geom_text(
        mapping       = aes(x = Inf, y = Inf, label = n),
        color         = "black",
        size          = lab_pt,
        alpha         = 1,
        check_overlap = TRUE,
        hjust         = 1.1,
        vjust         = v_jst
      )
  }
  
  res <- res +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, y_exp)),
      trans = y_trans
    )
  
  if (identical(n_trts, 1L)) {
    res <- res +
      theme(legend.position = "none")
  }
  
  res
}

create_stat_vlns <- function(df_in, metric, reg_labels = regs) {
  
  reg_labels <- reg_labels[names(reg_labels) %in% df_in$region]
  
  # Format plot data
  dat <- df_in %>%
    dplyr::filter(type == metric) %>%
    
    group_by(rep, type, region, sub_reg) %>%
    mutate(
      n     = n_distinct(name),
      n_lab = str_c("n = ", comma(n))
    ) %>%
    ungroup() %>%
    mutate(
      region = reg_labels[region],
      region = fct_relevel(region, unname(reg_labels))
    )
  
  # Create violins
  res <- dat %>%
    ggplot(aes(region, fc)) +
    
    geom_boxplot(notch = TRUE, alpha = 0.5, fill = last(clrs)) +    
    # geom_violin(draw_quantiles = c(0.25, 0.75), fill = last(clrs), alpha = 0.5) +
    # stat_summary(geom = "point", fun = median, size = 0.5) +
    geom_hline(yintercept = 1, linetype = 2) +
    
    geom_text(
      mapping       = aes(x = Inf, y = Inf, label = n_lab),
      color         = "black",
      size          = txt_pt / .pt,
      alpha         = 1,
      check_overlap = TRUE,
      hjust         = 1.1,
      vjust         = 1.4
    ) +
    
    labs(y = str_c("fold change\n", names(metric))) +
    scale_y_log10() +
    plt_x_scale +
    box_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  res
}

create_multistat_vlns <- function(df_in, metrics, reg_labels = regs) {
  
  reg_labels <- reg_labels[names(reg_labels) %in% df_in$region]
  
  dat <- df_in %>%
    dplyr::filter(type %in% names(metrics)) %>%
    dplyr::select(-all_of(treats)) %>%
    
    group_by(rep, type, region, sub_reg) %>%
    mutate(
      n     = n_distinct(name),
      n_lab = str_c("n = ", comma(n))
    ) %>%
    ungroup() %>%
    mutate(
      region = reg_labels[region],
      region = fct_relevel(region, unname(reg_labels)),
      type   = fct_relevel(type, names(metrics))
    )
  
  res <- dat %>%
    ggplot(aes(type, fc, fill = type)) +
    geom_violin(
      draw_quantiles = c(0.25, 0.75),
      fill           = last(clrs),
      alpha          = 0.5
    ) +
    
    geom_text(
      mapping       = aes(x = Inf, y = Inf, label = n_lab),
      color         = "black",
      size          = txt_pt / .pt,
      alpha         = 1,
      check_overlap = TRUE,
      hjust         = 1.1,
      vjust         = 1.4
    ) +
    
    stat_summary(geom = "point", fun = median, size = 0.5) +
    geom_hline(yintercept = 1, linetype = 2) +
    
    scale_y_log10() +
    
    scale_x_discrete(labels = metrics) +
    
    facet_wrap(~ rep, scales = "free_y") +
    labs(y = "fold change") +
    box_theme +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 0, hjust = c(0.6, 0.4))
    )
  
  res
}

create_stat_corr <- function(df_in, x, y, ttl = NULL,
                             plt_clrs = c("black", last(clrs)),
                             reg_labels = regs) {
  
  reg_labels <- reg_labels[names(reg_labels) %in% df_in$region]
  
  if (all(plt_clrs %in% c("black", "#000000"))) {
    plt_clrs[2] <- lighten(plt_clrs[2], 0.7)
  }
  
  # Check for NAs in data
  corr_dat <- df_in %>%
    pivot_wider(names_from = "type", values_from = "counts")
  
  if (!all(complete.cases(corr_dat[c(x, y)]))) {
    stop("NAs present in data.")
  }
  
  # Data for scatter plots
  # x-axis is log10 transformed
  corr_dat <- corr_dat %>%
    group_by(treat, rep, region) %>%
    mutate(
      n     = n_distinct(name),
      n     = str_c("n = ", comma(n)),
      corr  = list(broom::tidy(cor.test(log10(!!sym(x)), !!sym(y)), method = "spearman")),
      r     = map_dbl(corr, ~ .x$estimate),
      r     = round(r, 2),
      r     = str_c("r = ", r),
      n_lab = str_c(n, "\n", r)
    ) %>%
    ungroup() %>%
    mutate(
      treat  = fct_relevel(treat, treats),
      region = reg_labels[region],
      region = fct_relevel(region, unname(reg_labels))
    )
  
  # Create scatter plots
  res <- corr_dat %>%
    ggplot(aes(!!sym(x), !!sym(y))) +
    geom_pointdensity(size = 0.01) +
    geom_smooth(method = "lm", formula = y ~ x, color = "black", linetype = 2, size = 0.5) +
    
    geom_text(
      aes(Inf, Inf, label = n_lab),
      hjust = 1.07,
      vjust = 1.1,
      color = "black",
      size  = txt_pt / 1.2 / .pt,
      check_overlap = TRUE
    ) +
    
    facet_grid(rep ~ region) +
    
    scale_color_gradientn(colours = plt_clrs) +
    scale_x_log10() +
    guides(color = guide_colorbar(title = "n genes", barwidth = unit(5, "pt"))) +
    
    labs(title = ttl, x = names(x), y = names(y)) +
    theme_info +
    theme(
      aspect.ratio = 0.9,
      legend.title = element_text(size = ttl_pt),
      legend.text  = element_text(size = txt_pt),
      axis.line    = element_blank(),
      panel.border = element_rect(fill = NA, color = "black")
    )
  
  res
}

# Identify genes that have >1 pause in both exonic and intronic regions
get_ex_int_genes <- function(df_in, data_type, samples = unname(sam_nms), min_counts = 1,
                             len_norm = TRUE) {
  
  reg_filt <- unique(df_in$region)
  
  res <- df_in %>%
    dplyr::filter(type == data_type) %>%
    group_by(sample, region, gene_name) %>%
    summarize(
      len     = sum((end - start) / 1000),
      counts  = sum(counts),
      .groups = "drop"
    )
  
  if (len_norm) {
    res <- res %>%
      mutate(counts = counts / len)
  }
  
  res <- res %>%
    group_by(sample, gene_name) %>%
    dplyr::filter(
      all(counts > min_counts),
      all(reg_filt %in% region)
    ) %>%
    group_by(gene_name, region) %>%
    dplyr::filter(all(samples %in% sample)) %>%
    ungroup()
  
  stopifnot(n_distinct(table(res$region)) == 1)
  
  res <- unique(res$gene_name)
  
  res
}

# Sum intron and exon counts
sum_ex_int_counts <- function(df_in, sum_grps, n_grps, gene_list = NULL) {
  
  # Sum counts for all exon and intron regions for each gene
  # group_by(sample, treat, rep, type, gene_name, txn_id, region) %>%
  res <- df_in %>%
    group_by(!!!syms(sum_grps)) %>%
    summarize(
      counts  = sum(counts),
      len     = sum((end - start) / 1000),
      .groups = "drop"
    )
  
  # Calculate stats to plot
  # filter exons/introns based on NET-seq signal and only keep genes with counts
  # for both exons and introns
  stats     <- c("NET", "p-reads", "pauses")
  log_stats <- str_c(stats, "_kb")
  
  res <- res %>%
    pivot_wider(names_from = type, values_from = counts)
  
  if (!is.null(gene_list)) {
    res <- res %>%
      dplyr::filter(gene_name %in% gene_list)
  }
  
  res <- res %>%
    calc_pause_stats() %>%
    group_by(!!!syms(n_grps)) %>%
    mutate(
      n = n_distinct(gene_name),
      n = str_c("n = ", comma(n))
    ) %>%
    ungroup()
  
  res
}

# Plot pause stats for exons and introns
plot_ex_int_pauses <- function(df_in, y = "p-reads_NET", y_ttl = "fraction pause reads",
                               violin = TRUE, plot_clrs = clrs, p_lab = FALSE,
                               p_grps = "sample", ...) {
  
  dat <- df_in %>%
    dplyr::filter(type == y)
  
  # Add p-value labels
  if (p_lab) {
    lab_dat <- dat %>%
      dplyr::select(-len) %>%
      pivot_wider(names_from = "region", values_from = "counts") %>%
      group_by(!!!syms(c(p_grps, "n"))) %>%
      summarize(
        p = (wilcox.test(exons, introns))$p.value,
        p = ifelse(p == 0, "0", scientific(p, 2)),
        p = str_c(n, "\np = ", p),
        .groups = "drop"
      ) %>%
      distinct()
    
    dat <- dat %>%
      left_join(lab_dat, by = c(p_grps, "n")) %>%
      mutate(n = p)
  }
  
  # Create ggplot
  res <- dat %>%
    ggplot(aes(region, counts, alpha = region, fill = treat))
  
  if (violin) {
    res <- res +
      geom_violin(draw_quantiles = c(0.25, 0.75), scale = "width", ...) +
      stat_summary(geom = "point", fun = median, alpha = 1)
    
  } else {
    res <- res +
      geom_boxplot(notch = TRUE, outlier.size = 0.3, ...)
  }
  
  res <- res +    
    geom_text(
      aes(2.5, Inf, label = n),
      hjust = 1.1,
      vjust = 1.3,
      color = "black",
      alpha = 1,
      size  = txt_pt / .pt,
      check_overlap = TRUE
    ) +
    scale_alpha_manual(values = c(0.3, 0.6)) +
    scale_fill_manual(values = plot_clrs) +
    labs(y = y_ttl) +
    theme_info +
    theme(
      aspect.ratio    = 1.5,
      legend.position = "none",
      axis.title.x    = element_blank(),
      axis.text.x     = element_text(size = ttl_pt)
    )
  
  res
}
