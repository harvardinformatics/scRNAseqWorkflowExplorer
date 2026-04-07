jaccard_similarity <- function(set1, set2) {
  intersect_length <- length(intersect(set1, set2))
  union_length <- length(set1) + length(set2) - intersect_length
  intersect_length / union_length
}

MakeInterVsIntraStablePlot <- function(meta1, meta2,
                                       bootstraps1, bootstraps2,
                                       threshold, label1, label2) {

  clusters1 <- tibble::tibble(
    cellbarcode = rownames(meta1),
    clusterid1 = as.character(meta1$seurat_clusters)
  )

  clusters2 <- tibble::tibble(
    cellbarcode = rownames(meta2),
    clusterid2 = as.character(meta2$seurat_clusters)
  )

  clusters_merged <- dplyr::full_join(clusters1, clusters2, by = "cellbarcode")

  indices1 <- split(seq_len(nrow(clusters_merged)), clusters_merged$clusterid1)
  indices2 <- split(seq_len(nrow(clusters_merged)), clusters_merged$clusterid2)

  indices1 <- indices1[!is.na(names(indices1))]
  indices2 <- indices2[!is.na(names(indices2))]

  jaccard_list <- list()
  for (i in names(indices1)) {
    for (j in names(indices2)) {
      set1 <- indices1[[i]]
      set2 <- indices2[[j]]
      similarity <- jaccard_similarity(set1, set2)
      jaccard_list[[length(jaccard_list) + 1]] <- list(
        name1 = i,
        name2 = j,
        jaccard_similarity = similarity
      )
    }
  }

  jaccard_df <- do.call(rbind, lapply(jaccard_list, as.data.frame))
  jaccard_df <- type.convert(jaccard_df, as.is = TRUE)

  jaccard_df$name1 <- factor(jaccard_df$name1,
    levels = sort(unique(jaccard_df$name1))
  )
  jaccard_df$name2 <- factor(jaccard_df$name2,
    levels = sort(unique(jaccard_df$name2))
  )
  jaccard_tibble <- tibble::as_tibble(jaccard_df)

  stable_cluster_count <- jaccard_tibble %>%
    dplyr::filter(jaccard_similarity >= threshold) %>%
    dplyr::summarise(n = dplyr::n_distinct(name1)) %>%
    dplyr::pull(n)

  bootstrap1_summary <- bootstraps1 %>%
    dplyr::group_by(bootstrap_number) %>%
    dplyr::filter(max_jaccard >= threshold) %>%
    dplyr::summarise(n_clusters = dplyr::n_distinct(clusterid), .groups = "drop") %>%
    dplyr::mutate(method = "n_stable_name1")

  bootstrap2_summary <- bootstraps2 %>%
    dplyr::group_by(bootstrap_number) %>%
    dplyr::filter(max_jaccard >= threshold) %>%
    dplyr::summarise(n_clusters = dplyr::n_distinct(clusterid), .groups = "drop") %>%
    dplyr::mutate(method = "n_stable_name2")

  boot_stable_merged <- dplyr::bind_rows(bootstrap1_summary, bootstrap2_summary) %>%
    dplyr::mutate(method = factor(method,
      levels = c("n_stable_name1", "n_stable_name2")
    ))

  method_colors <- c(
    n_stable_name1 = "forestgreen",
    n_stable_name2 = "dodgerblue3"
  )

  max_y <- max(c(boot_stable_merged$n_clusters, stable_cluster_count), na.rm = TRUE)
  upper_y <- max(1, max_y) * 1.15
  y_break_step <- max(1, floor(max_y / 8))
  label_y <- -max(1, max_y) * 0.08
  label_data <- tibble::tibble(
    method = factor(c("n_stable_name1", "n_stable_name2"),
      levels = c("n_stable_name1", "n_stable_name2")
    ),
    x_pos = c(1.08, 1.92),
    y_pos = c(label_y * 0.85, label_y * 1.25),
    label = c(label1, label2)
  )

  boot_stable_merged %>%
    ggplot2::ggplot(ggplot2::aes(x = method, y = n_clusters, fill = method, color = method)) +
    ggplot2::geom_hline(yintercept = -0.5, color = "black", linewidth = 0.5) +
    ggplot2::geom_violin(width = 0.7, trim = FALSE, alpha = 0.35, linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = stable_cluster_count,
      color = "black", linetype = "dashed", linewidth = 1
    ) +
    ggplot2::annotate(
      "text",
      x = 1.5,
      y = stable_cluster_count,
      label = paste0(
        "number of inter-method stable clusters, min. Jaccard similarity threshold = ",
        as.character(threshold)
      ),
      vjust = -0.6,
      size = 5.2,
      color = "black"
    ) +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    ggplot2::geom_text(
      data = label_data,
      ggplot2::aes(x = x_pos, y = y_pos, label = label, color = method),
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = 1,
      size = 3.4
    ) +
    ggplot2::scale_x_discrete(
      labels = NULL,
      expand = c(0.20, 0.20)
    ) +
    ggplot2::scale_fill_manual(values = method_colors, guide = "none") +
    ggplot2::scale_color_manual(values = method_colors, guide = "none") +
    ggplot2::scale_y_continuous(
      breaks = seq(0, max_y, by = y_break_step),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::coord_cartesian(ylim = c(-0.5, upper_y), clip = "off") +
    ggplot2::labs(
      x = "",
      y = "# stable clusters"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = ggplot2::element_line(color = "black", linewidth = 0.5),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 46, 10)
    )
}

make_contrast_cluster_palette <- function(n_colors) {
  if (n_colors <= 0) {
    return(character())
  }

  if (n_colors == 1) {
    return("#1B9E77")
  }

  make_hcl_hex <- function(h, c, l) {
    grDevices::hcl(h = h %% 360, c = c, l = l)
  }

  anchor_hues <- c(15, 195)
  split_offsets <- c(-32, 32, -58, 58, -85, 85, -120, 120, -150, 150)
  luminance_levels <- c(52, 68)
  chroma_levels <- c(85, 72)

  candidate_colors <- c()
  for (lum in luminance_levels) {
    for (chr in chroma_levels) {
      for (anchor in anchor_hues) {
        candidate_colors <- c(candidate_colors, make_hcl_hex(anchor, chr, lum))
        for (offset in split_offsets) {
          candidate_colors <- c(candidate_colors, make_hcl_hex(anchor + offset, chr, lum))
        }
      }
    }
  }

  candidate_count <- max(24, n_colors * 8)
  candidate_colors <- unique(c(
    candidate_colors,
    grDevices::hcl.colors(candidate_count, palette = "Dynamic"),
    grDevices::hcl.colors(candidate_count, palette = "Dark 3")
  ))

  candidate_colors <- candidate_colors[seq_len(min(length(candidate_colors), max(candidate_count, n_colors)))]
  candidate_rgb <- t(grDevices::col2rgb(candidate_colors) / 255)
  candidate_lab <- grDevices::convertColor(candidate_rgb, from = "sRGB", to = "Lab", scale.in = 1)

  selected_idx <- integer(n_colors)
  selected_idx[[1]] <- 1L
  selected_idx[[2]] <- which.max(rowSums((candidate_lab - matrix(candidate_lab[selected_idx[[1]], ], nrow(candidate_lab), 3, byrow = TRUE))^2))

  available_idx <- setdiff(seq_len(nrow(candidate_lab)), selected_idx[seq_len(2)])

  if (n_colors > 2) {
    for (i in 3:n_colors) {
      recent_idx <- selected_idx[seq_len(i - 1)]
      recent_lab <- candidate_lab[recent_idx, , drop = FALSE]
      distance_matrix <- vapply(
        available_idx,
        function(idx) {
          rowSums((recent_lab - matrix(candidate_lab[idx, ], nrow(recent_lab), 3, byrow = TRUE))^2)
        },
        numeric(length(recent_idx))
      )

      if (is.null(dim(distance_matrix))) {
        distance_matrix <- matrix(distance_matrix, nrow = length(recent_idx))
      }

      candidate_scores <- apply(distance_matrix, 2, function(distances) {
        primary_gap <- distances[[length(distances)]]
        secondary_gap <- if (length(distances) >= 2) distances[[length(distances) - 1]] else primary_gap
        tertiary_gap <- if (length(distances) >= 3) distances[[length(distances) - 2]] else secondary_gap
        min(primary_gap, secondary_gap * 0.95, tertiary_gap * 0.8, min(distances) * 0.65)
      })

      next_idx <- available_idx[[which.max(candidate_scores)]]
      selected_idx[[i]] <- next_idx
      available_idx <- setdiff(available_idx, next_idx)
    }
  }

  candidate_colors[selected_idx]
}

make_expression_heatmap_plot <- function(plot_data, gene_symbol, cluster_annotations = tibble::tibble(), cluster_boundaries = numeric(), sort_mode = "expression", sort_method = NULL, for_pdf = FALSE) {
  method_levels <- levels(plot_data$method)
  n_methods <- length(method_levels)
  has_display_x <- "display_x" %in% names(plot_data)

  plot_data <- plot_data %>%
    dplyr::mutate(
      x_index = if (has_display_x) display_x else as.numeric(barcode),
      y_index = n_methods - as.integer(method) + 1
    )

  method_axis_labels <- if (!is.null(sort_method) && nzchar(sort_method) && sort_method %in% method_levels) {
    stats::setNames(
      ifelse(method_levels == sort_method, paste0(method_levels, "**"), method_levels),
      method_levels
    )
  } else {
    stats::setNames(method_levels, method_levels)
  }
  legend_data <- tibble::tibble(
    status_label = c("Barcode missing", "Expression missing/undefined"),
    x_index = 1,
    y_index = 1
  )

  y_axis_text_size <- if (for_pdf) 10 else 8
  plot_title_size <- if (for_pdf) 13 else 11
  tile_height <- if (for_pdf) 0.84 else 0.92
  tile_width <- 1
  top_margin <- 48
  bottom_margin <- if (identical(sort_mode, "cluster") && nrow(cluster_annotations) > 0) 84 else 42
  boundary_bands <- if (identical(sort_mode, "cluster") && length(cluster_boundaries) > 0) {
    tibble::tibble(x = cluster_boundaries)
  } else {
    tibble::tibble()
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_index, y = y_index, fill = expression)) +
    ggplot2::geom_tile(width = tile_width, height = tile_height) +
    ggplot2::scale_fill_gradientn(
      colours = c("#313695", "#4575b4", "lightgoldenrod1", "#fdae61", "#f46d43", "#d73027"),
      na.value = "transparent",
      name = paste(gene_symbol, "Normalized gene expression", sep = "\n")
    ) +
    ggplot2::geom_tile(
      data = dplyr::filter(plot_data, cell_status == "barcode_missing"),
      fill = "black",
      inherit.aes = FALSE,
      ggplot2::aes(x = x_index, y = y_index),
      width = tile_width,
      height = tile_height
    ) +
    ggplot2::geom_tile(
      data = dplyr::filter(plot_data, cell_status == "expression_missing"),
      fill = "gray75",
      inherit.aes = FALSE,
      ggplot2::aes(x = x_index, y = y_index),
      width = tile_width,
      height = tile_height
    ) +
    ggplot2::geom_point(
      data = legend_data,
      ggplot2::aes(x = x_index, y = y_index, color = status_label),
      inherit.aes = FALSE,
      alpha = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Barcode missing" = "black",
        "Expression missing/undefined" = "gray75"
      ),
      name = "Cell status"
    ) +
    ggplot2::labs(
      title = NULL,
      x = "",
      y = NULL
    ) +
    ggplot2::scale_x_continuous(
      breaks = NULL,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(method_levels),
      labels = rev(unname(method_axis_labels)),
      expand = ggplot2::expansion(mult = c(if (identical(sort_mode, "cluster") && nrow(cluster_annotations) > 0) 0.24 else 0.02, 0.18))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      axis.text.y = ggplot2::element_text(size = y_axis_text_size, color = "black"),
      plot.title = ggplot2::element_text(size = plot_title_size),
      legend.title = ggplot2::element_text(size = if (for_pdf) 10 else 8, face = "bold"),
      legend.text = ggplot2::element_text(size = if (for_pdf) 9 else 7),
      panel.grid = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(fill = "white", color = "gray85"),
      plot.margin = ggplot2::margin(top_margin, 10, bottom_margin, 10)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        order = 2,
        override.aes = list(alpha = 1, shape = 22, size = 5, fill = c("black", "gray75"))
      ),
      fill = ggplot2::guide_colorbar(order = 1)
    ) +
    ggplot2::annotate(
      "text",
      x = mean(range(plot_data$x_index)),
      y = n_methods + 0.88,
      label = "Cell barcode",
      size = if (for_pdf) 4 else 3
    ) +
    ggplot2::annotate(
      "text",
      x = mean(range(plot_data$x_index)),
      y = if (identical(sort_mode, "cluster") && nrow(cluster_annotations) > 0) -1.06 else -0.5,
      label = "** method used to sort expression",
      size = if (for_pdf) 3.6 else 2.8
    )

  if (identical(sort_mode, "cluster") && nrow(cluster_annotations) > 0) {
    p <- p +
      ggplot2::geom_rect(
        data = cluster_annotations,
        ggplot2::aes(
          xmin = x_start - tile_width / 2,
          xmax = x_end + tile_width / 2,
          ymin = -0.08,
          ymax = 0.36
        ),
        inherit.aes = FALSE,
        fill = cluster_annotations$cluster_color,
        color = NA,
        show.legend = FALSE
      ) +
      ggplot2::geom_segment(
        data = boundary_bands,
        ggplot2::aes(
          x = x,
          xend = x,
          y = -0.08,
          yend = n_methods + tile_height / 2
        ),
        inherit.aes = FALSE,
        color = "white",
        linewidth = if (for_pdf) 1.5 else 2.2,
        lineend = "butt",
        show.legend = FALSE
      ) +
      ggplot2::annotate(
        "text",
        x = min(plot_data$x_index) - max(10, 0.015 * max(plot_data$x_index)),
        y = 0.14,
        label = "Cluster",
        hjust = 1,
        size = if (for_pdf) 4 else 3
      )
  }

  p
}

make_cluster_preview_plot <- function(plot_data, cluster_info, gene_symbol, sort_method = NULL) {
  method_levels <- levels(plot_data$method)
  has_display_x <- "display_x" %in% names(plot_data)
  preview_data <- plot_data %>%
    dplyr::mutate(
      x_index = if (has_display_x) display_x else as.numeric(barcode)
    ) %>%
    dplyr::filter(.data$x_index >= cluster_info$xmin, .data$x_index <= cluster_info$xmax) %>%
    dplyr::mutate(barcode_index = dplyr::dense_rank(.data$x_index))

  if (nrow(preview_data) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::annotate("text", x = 0, y = 0, label = "No cells in cluster")
    )
  }

  preview_data <- preview_data %>%
    dplyr::mutate(
      y_index = length(method_levels) - as.integer(method) + 1
    )

  ggplot2::ggplot(preview_data, ggplot2::aes(x = barcode_index, y = y_index, fill = expression)) +
    ggplot2::geom_tile(width = 1, height = 0.86) +
    ggplot2::geom_tile(
      data = dplyr::filter(preview_data, cell_status == "barcode_missing"),
      fill = "black",
      inherit.aes = FALSE,
      ggplot2::aes(x = barcode_index, y = y_index),
      width = 1,
      height = 0.86
    ) +
    ggplot2::geom_tile(
      data = dplyr::filter(preview_data, cell_status == "expression_missing"),
      fill = "gray75",
      inherit.aes = FALSE,
      ggplot2::aes(x = barcode_index, y = y_index),
      width = 1,
      height = 0.86
    ) +
    ggplot2::scale_fill_gradientn(
      colours = c("#313695", "#4575b4", "lightgoldenrod1", "#fdae61", "#f46d43", "#d73027"),
      na.value = "transparent",
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      breaks = NULL,
      labels = NULL,
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      title = paste(gene_symbol, "| cluster", cluster_info$label[[1]]),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 9, face = "bold"),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
}
