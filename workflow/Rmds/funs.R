# scripts from Ryan

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
      axis.text.x     = element_text(face ="bold", angle = 45, hjust = 1)
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
        face          ="bold",
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

