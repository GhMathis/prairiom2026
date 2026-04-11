# ============================================================
# plot_multipartite_sbm()
# Plot rearranged interaction matrices from a MultipartiteSBM_fit
# object, with observations ordered by their estimated cluster.
#
# Dependencies: ggplot2, dplyr, tidyr, patchwork
# ============================================================

plot_multipartite_sbm <- function(msbm) {
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(patchwork)
  
  # ----------------------------------------------------------
  # 1.  Memberships and per-FG sorted indices
  # ----------------------------------------------------------
  memberships <- msbm$memberships
  
  sorted_idx <- lapply(memberships, order)
  sorted_mem <- lapply(names(memberships), function(fg)
    memberships[[fg]][sorted_idx[[fg]]])
  names(sorted_mem) <- names(memberships)
  
  # ----------------------------------------------------------
  # 2.  Helper: cluster-separator positions along one axis
  # ----------------------------------------------------------
  separators <- function(mem) which(diff(mem) != 0) + 0.5
  
  # ----------------------------------------------------------
  # 3.  Helper: build one ggplot tile-matrix for one network
  #     is_last = TRUE  -> show y-axis labels and title
  #     is_last = FALSE -> hide them
  # ----------------------------------------------------------
  one_plot <- function(net, is_last = FALSE) {
    
    dl     <- net$dimLabels
    mat    <- net$networkData
    row_fg <- dl[1]
    col_fg <- dl[2]
    
    r_ord  <- sorted_idx[[row_fg]]
    c_ord  <- sorted_idx[[col_fg]]
    mat_r  <- mat[r_ord, c_ord, drop = FALSE]
    
    nr <- nrow(mat_r);  nc <- ncol(mat_r)
    
    col_nms <- if (!is.null(colnames(mat_r))) colnames(mat_r) else paste0("C", seq_len(nc))
    row_nms <- if (!is.null(rownames(mat_r))) rownames(mat_r) else paste0("R", seq_len(nr))
    
    df <- expand.grid(row_pos = seq_len(nr), col_pos = seq_len(nc)) %>%
      mutate(
        value = as.vector(mat_r),
        y     = nr - row_pos + 1
      )
    
    row_sep   <- separators(sorted_mem[[row_fg]])
    col_sep   <- separators(sorted_mem[[col_fg]])
    row_sep_y <- nr - row_sep + 1
    
    is_binary <- all(mat_r %in% c(0L, 1L, NA))
    
    fill_scale <- if (is_binary) {
      scale_fill_gradient(low = "#f7fbff", high = "#08306b",
                          na.value = "grey85", guide = "none")
    } else {
      scale_fill_gradientn(
        colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
        na.value = "grey85", name = "")
    }
    
    p <- ggplot(df, aes(x = col_pos, y = y, fill = value)) +
      geom_tile(linewidth = 0.08) +
      fill_scale +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(
        breaks   = if (is_last) seq_len(nr) else NULL,
        labels   = if (is_last) row_nms[nr:1] else NULL,
        expand   = c(0, 0),
        position = "right"
      ) +
      geom_hline(yintercept = row_sep_y,
                 colour = "#e31a1c", linewidth = 0.65) +
      geom_vline(xintercept = col_sep,
                 colour = "#e31a1c", linetype = "dashed", linewidth = 0.65) +
      labs(
        title = if (is_last) paste0(col_fg) else NULL,
        x     = col_fg,
        y     = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y  = if (is_last)
          element_text(face = "italic", colour = "grey30", size = 7)
        else
          element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid   = element_blank(),
        panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.6),
        plot.title   = if (is_last)
          element_text(face = "bold", hjust = 0.5, size = 10, margin = margin(b = 4))
        else
          element_blank(),
        axis.title.x = element_text(face = "italic", colour = "grey30", size = 9),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width  = unit(0.25, "cm"),
        legend.text       = element_text(size = 7)
      )
    
    p
  }
  
  # ----------------------------------------------------------
  # 4.  Build one plot per network — flag the last one
  # ----------------------------------------------------------
  n_nets <- length(msbm$networkData)
  
  plots <- lapply(seq_len(n_nets), function(i)
    one_plot(msbm$networkData[[i]], is_last = (i == n_nets)))
  
  wrap_plots(plots, nrow = 1) +
    plot_annotation(
      title   = "Multipartite SBM — matrices ordered by block",
      caption = "Red lines delimit estimated blocks",
      theme   = theme(
        plot.title   = element_text(face = "bold", size = 13, hjust = 0.5,
                                    margin = margin(b = 6)),
        plot.caption = element_text(colour = "#e31a1c", size = 8, hjust = 0.5)
      )
    )
}


plot_bipartite_sbm <- function(bsbm = sbm_Bact_bin, mat = NULL) {
  
  require(ggplot2)
  require(dplyr)
  require(patchwork)
  
  # ----------------------------------------------------------
  # 1.  Labels and memberships
  # ----------------------------------------------------------
  dl      <- bsbm$dimLabels          # c(row = "Grids", col = "Plants")
  row_fg  <- dl[1]
  col_fg  <- dl[2]
  
  mem_row <- bsbm$memberships[[row_fg]]
  mem_col <- bsbm$memberships[[col_fg]]
  
  # Sort indices by cluster label
  r_ord   <- order(mem_row)
  c_ord   <- order(mem_col)
  
  sorted_mem_row <- mem_row[r_ord]
  sorted_mem_col <- mem_col[c_ord]
  
  # ----------------------------------------------------------
  # 2.  Reorder the matrix
  # ----------------------------------------------------------
  
  if(is.null(mat)){
    mat <- bsbm$networkData
  }
  mat_r <- mat[r_ord, c_ord, drop = FALSE]
  
  nr <- nrow(mat_r)
  nc <- ncol(mat_r)
  
  row_nms <- if (!is.null(rownames(mat_r))) rownames(mat_r) else paste0("R", seq_len(nr))
  col_nms <- if (!is.null(colnames(mat_r))) colnames(mat_r) else paste0("C", seq_len(nc))
  
  # ----------------------------------------------------------
  # 3.  Long format for ggplot
  # ----------------------------------------------------------
  df <- expand.grid(row_pos = seq_len(nr), col_pos = seq_len(nc)) %>%
    mutate(
      value = as.vector(mat_r),
      y     = nr - row_pos + 1          # flip: row 1 at top
    )
  
  # ----------------------------------------------------------
  # 4.  Cluster separator positions
  # ----------------------------------------------------------
  separators <- function(mem) which(diff(mem) != 0) + 0.5
  
  row_sep   <- separators(sorted_mem_row)
  col_sep   <- separators(sorted_mem_col)
  row_sep_y <- nr - row_sep + 1         # convert to flipped-y space
  
  # ----------------------------------------------------------
  # 5.  Colour scale
  # ----------------------------------------------------------
  is_binary <- all(mat_r %in% c(0L, 1L, NA))
  
  fill_scale <- if (is_binary) {
    scale_fill_gradient(low  = "#f7fbff", high = "#08306b",
                        na.value = "grey85", guide = "none")
  } else {
    scale_fill_gradientn(
      colours  = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
      na.value = "grey85",
      name     = "")
  }
  
  # ----------------------------------------------------------
  # 6.  Block-size annotations for axis strips
  # ----------------------------------------------------------
  row_block_sizes <- table(sorted_mem_row)
  col_block_sizes <- table(sorted_mem_col)
  
  # ----------------------------------------------------------
  # 7.  Plot
  # ----------------------------------------------------------
  p <- ggplot(df, aes(x = col_pos, y = y, fill = value)) +
    geom_tile( linewidth = 0.08) +
    fill_scale +
    # x-axis: hidden (columns)
    scale_x_continuous(expand = c(0, 0)) +
    # y-axis: row names on the right
    scale_y_continuous(
      breaks   = seq_len(nr),
      labels   = row_nms[nr:1],
      expand   = c(0, 0),
      position = "right"
    ) +
    # Block separator lines
    geom_hline(yintercept = row_sep_y,
               colour = "#e31a1c", linewidth = 0.7) +
    geom_vline(xintercept = col_sep,
               colour = "#e31a1c", linetype = "dashed", linewidth = 0.7) +
    labs(
      title    = paste0(row_fg, "  ×  ", col_fg),
      subtitle = paste0(
        length(row_block_sizes), " ", row_fg, " blocks  |  ",
        length(col_block_sizes), " ", col_fg, " blocks"
      ),
      x = col_fg,
      y = row_fg
    ) +
    theme_minimal(base_size = 11) +
    theme(
      # x-axis (columns): hidden — too many
      axis.text.x   = element_blank(),
      axis.ticks.x  = element_blank(),
      # y-axis (rows): shown, italic, small
      axis.text.y   = element_text(face = "italic", colour = "grey30", size = 7),
      axis.ticks.y  = element_blank(),
      # clean background
      panel.grid    = element_blank(),
      panel.border  = element_rect(colour = "grey40", fill = NA, linewidth = 0.7),
      # titles
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13,
                                   margin = margin(b = 3)),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40",
                                   margin = margin(b = 6)),
      # axis labels
      axis.title.x  = element_text(face = "italic", colour = "grey30", size = 10),
      axis.title.y  = element_text(face = "italic", colour = "grey30", size = 10),
      # legend
      legend.position   = "left",
      legend.key.height = unit(1, "cm"),
      legend.key.width  = unit(0.3, "cm"),
      legend.text       = element_text(size = 8)
    )
  
  # Wrap with caption
  p + plot_annotation(
    caption = "Red lines delimit estimated blocks",
    theme   = theme(
      plot.caption = element_text(colour = "#e31a1c", size = 8, hjust = 0.5)
    )
  )
}

# ============================================================
# plot_bipartite_sbm_avg()
# Compute the average SBM expectation per cluster along one
# axis, then plot the resulting compressed matrix.
#
# by = "row" → keep individual rows, average cols into col-clusters
# by = "col" → keep individual cols, average rows into row-clusters
#
# Returns a list: $plot and $matrix
#
# Dependencies: ggplot2, dplyr, tidyr, tibble, patchwork
# ============================================================

plot_bipartite_sbm_avg <- function(bsbm, by = "col") {
  
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(tibble)
  require(patchwork)
  
  # ----------------------------------------------------------
  # 1.  Labels and memberships
  # ----------------------------------------------------------
  dl     <- bsbm$dimLabels
  row_fg <- dl[1]   # e.g. "grid"
  col_fg <- dl[2]   # e.g. "bact"
  
  mem_row <- bsbm$memberships[[row_fg]]   # int vector, length = nrow
  mem_col <- bsbm$memberships[[col_fg]]   # int vector, length = ncol
  
  net     <- bsbm$networkData             # raw matrix (rows × cols)
  exp_mat <- bsbm$expectation             # same dims, fitted values
  
  row_nms_all <- rownames(net)
  col_nms_all <- colnames(net)
  
  # ----------------------------------------------------------
  # 2.  Build averaged matrix
  #     by = "col"  → rows stay individual, cols averaged per cluster
  #     by = "row"  → cols stay individual, rows averaged per cluster
  # ----------------------------------------------------------
  avg_mat <- if (by == "col") {
    
    # Average over col-clusters, keep rows as-is
    as.data.frame(exp_mat) %>%
      rownames_to_column(row_fg) %>%
      pivot_longer(-all_of(row_fg),
                   names_to  = col_fg,
                   values_to = "value") %>%
      left_join(setNames(data.frame(mem_col, col_nms_all), c("cluster", col_fg)),
                by = col_fg) %>%
      group_by(across(all_of(row_fg)), cluster) %>%
      summarise(value = mean(value), .groups = "drop") %>%
      mutate(cluster = paste0(col_fg, "_", cluster)) %>%
      pivot_wider(names_from = "cluster", values_from = "value") %>%
      column_to_rownames(row_fg) %>%
      as.matrix()
    
  } else {
    
    # Average over row-clusters, keep cols as-is
    as.data.frame(t(exp_mat)) %>%
      rownames_to_column(col_fg) %>%
      pivot_longer(-all_of(col_fg),
                   names_to  = row_fg,
                   values_to = "value") %>%
      left_join(setNames(data.frame(mem_row, row_nms_all), c("cluster", row_fg)),
                by = row_fg) %>%
      group_by(across(all_of(col_fg)), cluster) %>%
      summarise(value = mean(value), .groups = "drop") %>%
      mutate(cluster = paste0(row_fg, "_", cluster)) %>%
      pivot_wider(names_from = "cluster", values_from = "value") %>%
      column_to_rownames(col_fg) %>%
      as.matrix()
  }
  
  # ----------------------------------------------------------
  # 3.  Reorder averaged matrix by cluster memberships
  # ----------------------------------------------------------
  if (by == "col") {
    # rows ordered by row-cluster, cols already named by cluster number
    r_ord       <- order(mem_row)
    sorted_rows <- mem_row[r_ord]
    avg_mat     <- avg_mat[r_ord, , drop = FALSE]
    
    # Order col-cluster columns numerically
    col_cluster_nums <- as.integer(gsub(paste0(col_fg, "_"), "", colnames(avg_mat)))
    c_ord       <- order(col_cluster_nums)
    avg_mat     <- avg_mat[, c_ord, drop = FALSE]
    sorted_cols <- sort(col_cluster_nums)
    
  } else {
    # cols ordered by col-cluster, rows already named by cluster number
    c_ord       <- order(mem_col)
    sorted_cols <- mem_col[c_ord]
    avg_mat     <- avg_mat[, c_ord, drop = FALSE]
    
    row_cluster_nums <- as.integer(gsub(paste0(row_fg, "_"), "", rownames(avg_mat)))
    r_ord       <- order(row_cluster_nums)
    avg_mat     <- avg_mat[r_ord, , drop = FALSE]
    sorted_rows <- sort(row_cluster_nums)
  }
  
  nr <- nrow(avg_mat)
  nc <- ncol(avg_mat)
  
  row_nms <- rownames(avg_mat)
  col_nms <- colnames(avg_mat)
  
  # ----------------------------------------------------------
  # 4.  Separator helper and positions
  # ----------------------------------------------------------
  separators <- function(mem) which(diff(mem) != 0) + 0.5
  
  row_sep   <- separators(sorted_rows)
  col_sep   <- separators(sorted_cols)
  row_sep_y <- nr - row_sep + 1
  
  # ----------------------------------------------------------
  # 5.  Long-format data frame
  # ----------------------------------------------------------
  df <- expand.grid(row_pos = seq_len(nr), col_pos = seq_len(nc)) %>%
    mutate(
      value = as.vector(avg_mat),
      y     = nr - row_pos + 1
    )
  
  # ----------------------------------------------------------
  # 6.  Colour scale (always continuous for expectations)
  # ----------------------------------------------------------
  fill_scale <- scale_fill_gradientn(
    colours  = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    na.value = "grey85",
    name     = "Mean\nexpectation"
  )
  
  # ----------------------------------------------------------
  # 7.  Axis labels: shown on both axes (matrix is compressed)
  # ----------------------------------------------------------
  show_row_labels <- (by == "col")   # rows are real samples → show them
  show_col_labels <- (by == "row")   # cols are real samples → show them
  
  x_labels <- if (show_col_labels) col_nms else
    gsub(paste0(col_fg, "_"), "Cl. ", col_nms)   # cluster numbers on x
  
  y_labels <- if (show_row_labels) row_nms[nr:1] else
    gsub(paste0(row_fg, "_"), "Cl. ", rev(row_nms))
  
  # ----------------------------------------------------------
  # 8.  Number of blocks in each dimension
  # ----------------------------------------------------------
  n_row_blocks <- length(unique(mem_row))
  n_col_blocks <- length(unique(mem_col))
  
  # ----------------------------------------------------------
  # 9.  Plot
  # ----------------------------------------------------------
  p <- ggplot(df, aes(x = col_pos, y = y, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.15) +
    fill_scale +
    scale_x_continuous(
      breaks = seq_len(nc),
      labels = x_labels,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks   = seq_len(nr),
      labels   = y_labels,
      expand   = c(0, 0),
      position = "right"
    ) +
    geom_hline(yintercept = row_sep_y,
               colour = "#e31a1c", linewidth = 0.7) +
    geom_vline(xintercept = col_sep,
               colour = "#e31a1c", linetype = "dashed", linewidth = 0.7) +
    labs(
      title    = paste0(row_fg, "  ×  ", col_fg, "  (avg. expectation by ",
                        ifelse(by == "col", col_fg, row_fg), " cluster)"),
      subtitle = paste0(n_row_blocks, " ", row_fg, " blocks  |  ",
                        n_col_blocks, " ", col_fg, " blocks"),
      x = col_fg,
      y = row_fg
    ) +
    theme_minimal(base_size = 11) +
    theme(
      # x-axis: small rotated labels when showing cluster numbers,
      #         hidden when showing individual samples
      axis.text.x  = if (show_col_labels)
        element_text(angle = 90, vjust = 0.5, hjust = 1,
                     face = "italic", colour = "grey30", size = 6)
      else
        element_text(face = "bold", colour = "grey30", size = 9),
      axis.ticks.x = element_blank(),
      
      # y-axis: small italic when showing individual samples,
      #         bold cluster numbers otherwise
      axis.text.y  = if (show_row_labels)
        element_text(face = "italic", colour = "grey30", size = 7)
      else
        element_text(face = "bold", colour = "grey30", size = 9),
      axis.ticks.y = element_blank(),
      
      panel.grid   = element_blank(),
      panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.7),
      
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 12,
                                   margin = margin(b = 3)),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40",
                                   margin = margin(b = 6)),
      
      axis.title.x = element_text(face = "italic", colour = "grey30", size = 10),
      axis.title.y = element_text(face = "italic", colour = "grey30", size = 10),
      
      legend.position   = "left",
      legend.key.height = unit(1.2, "cm"),
      legend.key.width  = unit(0.35, "cm"),
      legend.text       = element_text(size = 8),
      legend.title      = element_text(size = 9, face = "bold")
    )
  
  final_plot <- p + plot_annotation(
    caption = "Red lines delimit estimated blocks",
    theme   = theme(
      plot.caption = element_text(colour = "#e31a1c", size = 8, hjust = 0.5)
    )
  )
  
  # ----------------------------------------------------------
  # 10.  Return both the plot and the averaged matrix
  # ----------------------------------------------------------
  invisible(list(
    plot   = final_plot,
    matrix = avg_mat
  ))
}