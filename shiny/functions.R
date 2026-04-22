jaccard_similarity <- function(set1, set2) {
  intersect_length <- length(intersect(set1, set2))
  union_length <- length(set1) + length(set2) - intersect_length
  intersect_length / union_length
}

make_shared_barcodes_upset_plot <- function(seurat_objects, method_names = NULL, min_size = 1, for_pdf = FALSE) {
  if (is.null(names(seurat_objects)) && is.null(method_names)) {
    stop("Provide `method_names` or a named list of Seurat objects.")
  }

  if (!is.list(seurat_objects) || length(seurat_objects) == 0) {
    stop("`seurat_objects` must be a non-empty list of Seurat objects.")
  }

  if (is.null(method_names)) {
    method_names <- names(seurat_objects)
  }

  if (length(method_names) != length(seurat_objects)) {
    stop("`method_names` must have the same length as `seurat_objects`.")
  }

  if (length(unique(method_names)) != length(method_names)) {
    stop("`method_names` must be unique.")
  }

  if (!is.numeric(min_size) || length(min_size) != 1 || is.na(min_size) || min_size < 1) {
    stop("`min_size` must be a single number greater than or equal to 1.")
  }

  min_size <- as.integer(min_size)

  invalid_objects <- !purrr::map_lgl(seurat_objects, inherits, "Seurat")
  if (any(invalid_objects)) {
    stop(
      "All entries in `seurat_objects` must inherit from Seurat. Invalid methods: ",
      paste(method_names[invalid_objects], collapse = ", ")
    )
  }

  barcode_membership <- purrr::map2_dfr(seurat_objects, method_names, function(obj, method_name) {
    tibble::tibble(
      cell_barcode = colnames(obj),
      method = method_name
    )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::distinct(cell_barcode, method) %>%
    dplyr::mutate(member = TRUE) %>%
    tidyr::pivot_wider(
      names_from = method,
      values_from = member,
      values_fill = FALSE
    )

  max_intersections <- min(max(1, 2^length(method_names) - 1), 20)
  base_text_size <- if (for_pdf) 11 else 12
  intersection_sizes <- barcode_membership %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      intersection_size = sum(dplyr::c_across(dplyr::all_of(method_names)))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(intersection_size > 0) %>%
    dplyr::count(dplyr::across(dplyr::all_of(method_names)), name = "cells", .drop = FALSE) %>%
    dplyr::filter(cells >= min_size) %>%
    dplyr::arrange(dplyr::desc(cells))

  annotation_y <- if (nrow(intersection_sizes) > 0) {
    max(intersection_sizes$cells) * 0.98
  } else {
    1
  }
  ComplexUpset::upset(
    barcode_membership,
    intersect = method_names,
    name = "Methods",
    sort_sets = FALSE,
    themes = list(
      overall = ggplot2::theme_bw(base_size = base_text_size) +
        ggplot2::theme(
          axis.title.x = ggplot2::element_text(),
          plot.title = ggplot2::element_text(face = "bold"),
          legend.position = "right",
          panel.grid.minor = ggplot2::element_blank()
        ),
      intersections_matrix = ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank()
      ),
      "Shared barcodes" = ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank()
      ),
      overall_sizes = ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
        axis.ticks.y = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(l = 18)
      )
    ),
    width_ratio = 0.18,
    min_size = min_size,
    n_intersections = max_intersections,
    base_annotations = list(
      "Shared barcodes" = ComplexUpset::intersection_size(
        counts = FALSE
      ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::annotate(
          geom = "text",
          x = Inf,
          y = annotation_y,
          label = paste0("minimum intersection size = ", scales::comma(min_size)),
          color = "black",
          hjust = 1.12,
          vjust = 1,
          size = if (for_pdf) 3.2 else 3
        ) +
        ggplot2::annotation_custom(
          grob = grid::textGrob(
            "# of cell barcodes",
            x = grid::unit(if (for_pdf) -3.2 else -2.9, "lines"),
            y = grid::unit(0.5, "npc"),
            rot = 90,
            just = "center",
            gp = grid::gpar(fontsize = if (for_pdf) 11 else 10)
          ),
          xmin = -Inf,
          xmax = -Inf,
          ymin = -Inf,
          ymax = Inf
        ) +
        ggplot2::labs(y = NULL)
    ),
    set_sizes = ComplexUpset::upset_set_size() +
      ggplot2::labs(y = "Method total")
  ) +
    ggplot2::labs(
      title = "Shared cell barcodes across methods",
      subtitle = paste0(length(method_names), " methods; ", scales::comma(nrow(barcode_membership)), " unique barcodes"),
      x = NULL
    )
}

arrange_cluster_levels <- function(cluster_values) {
  unique_clusters <- unique(cluster_values)
  suppressWarnings(cluster_numeric <- as.numeric(unique_clusters))

  if (all(!is.na(cluster_numeric))) {
    unique_clusters[order(cluster_numeric)]
  } else {
    sort(unique_clusters)
  }
}

jaccard_heatmap_plot <- function(meta1, meta2,
                                 name1, name2, threshold = 0.6,
                                 for_pdf = FALSE) {
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
      similarity <- jaccard_similarity(indices1[[i]], indices2[[j]])
      jaccard_list[[length(jaccard_list) + 1]] <- list(
        cluster1 = i,
        cluster2 = j,
        jaccard_similarity = similarity
      )
    }
  }

  jaccard_df <- do.call(rbind, lapply(jaccard_list, as.data.frame))
  jaccard_df <- type.convert(jaccard_df, as.is = TRUE)

  jaccard_df$cluster1 <- factor(
    jaccard_df$cluster1,
    levels = arrange_cluster_levels(as.character(jaccard_df$cluster1))
  )
  jaccard_df$cluster2 <- factor(
    jaccard_df$cluster2,
    levels = arrange_cluster_levels(as.character(jaccard_df$cluster2))
  )

  label_df <- jaccard_df %>%
    dplyr::filter(jaccard_similarity >= threshold) %>%
    dplyr::mutate(
      jaccard_label = if (for_pdf) "*" else format(round(jaccard_similarity, 2), nsmall = 2)
    )

  threshold_note <- paste0(
    "Cluster Jaccard stability\n",
    "threshold: ",
    format(threshold, trim = TRUE)
  )
  legend_note_df <- tibble::tibble(
    cluster1 = levels(jaccard_df$cluster1)[[1]],
    cluster2 = levels(jaccard_df$cluster2)[[1]],
    threshold_note = threshold_note
  )

  ggplot2::ggplot(
    data = jaccard_df,
    ggplot2::aes(x = cluster1, y = cluster2, fill = jaccard_similarity)
  ) +
    ggplot2::geom_tile(color = "black", linewidth = 0.35) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "firebrick",
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1)
    ) +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(label = jaccard_label),
      size = if (for_pdf) 4.2 else 2.8,
      color = "black",
      hjust = 0.5,
      vjust = if (for_pdf) 0.62 else 0.5
    ) +
    ggplot2::geom_point(
      data = legend_note_df,
      ggplot2::aes(x = cluster1, y = cluster2, alpha = threshold_note),
      inherit.aes = FALSE,
      shape = 15,
      size = 0,
      show.legend = TRUE
    ) +
    ggplot2::scale_x_discrete(name = name1) +
    ggplot2::scale_y_discrete(name = name2) +
    ggplot2::scale_alpha_manual(
      values = stats::setNames(0, threshold_note),
      guide = ggplot2::guide_legend(
        order = 2,
        title = NULL,
        override.aes = list(alpha = 0, size = 0)
      )
    ) +
    ggplot2::labs(fill = "Jaccard similarity") +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = if (for_pdf) 9 else 10, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = ggplot2::element_text(size = if (for_pdf) 9 else 10),
      axis.title.x = ggplot2::element_text(size = if (for_pdf) 11.2 else 14),
      axis.title.y = ggplot2::element_text(size = if (for_pdf) 11.2 else 14),
      legend.title.align = 0,
      legend.text.align = 0,
      legend.box = "vertical",
      legend.spacing.y = grid::unit(6, "pt")
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(order = 1),
      alpha = ggplot2::guide_legend(
        order = 2,
        title = NULL,
        label.hjust = 0,
        keywidth = grid::unit(0, "pt"),
        keyheight = grid::unit(0, "pt"),
        label.theme = ggplot2::element_text(hjust = 0, margin = ggplot2::margin(l = -8)),
        default.unit = "pt",
        override.aes = list(alpha = 0, size = 0)
      )
    )
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

  max_y <- max(c(boot_stable_merged$n_clusters, stable_cluster_count), na.rm = TRUE)
  annotation_offset <- max(0.75, max_y * 0.10)
  annotation_y <- max_y + annotation_offset

  p_value_data <- boot_stable_merged %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      p_value = {
        test_values <- n_clusters[!is.na(n_clusters)]
        if (length(test_values) == 0) {
          NA_real_
        } else {
          tryCatch(
            stats::wilcox.test(test_values, mu = stable_cluster_count, exact = FALSE)$p.value,
            error = function(e) NA_real_
          )
        }
      },
      y_pos = max(n_clusters, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(p_value), p_value <= 0.05) %>%
    dplyr::mutate(
      y_pos = annotation_y,
      label = paste0("p=", format(signif(p_value, 3), scientific = TRUE, trim = TRUE))
    )

  method_colors <- c(
    n_stable_name1 = "forestgreen",
    n_stable_name2 = "dodgerblue3"
  )

  upper_y <- max(
    max(1, max_y) + (2 * annotation_offset),
    if (nrow(p_value_data) > 0) max(p_value_data$y_pos) * 1.05 else 0
  )
  y_break_step <- max(1, floor(max_y / 8))
  inter_method_legend_label <- paste0(
    "Inter-method stable clusters\n(Jaccard threshold = ",
    format(threshold, trim = TRUE),
    ")"
  )
  inter_method_line_data <- tibble::tibble(
    yintercept = stable_cluster_count,
    legend_label = inter_method_legend_label
  )
  label_y <- -max(1, max_y) * 0.08
  label_data <- tibble::tibble(
    method = factor(c("n_stable_name1", "n_stable_name2"),
      levels = c("n_stable_name1", "n_stable_name2")
    ),
    x_pos = c(1.08, 1.92),
    y_pos = c(label_y, label_y),
    label = c(label1, label2)
  )

  boot_stable_merged %>%
    ggplot2::ggplot(ggplot2::aes(x = method, y = n_clusters, fill = method, color = method)) +
    ggplot2::geom_hline(yintercept = -0.5, color = "black", linewidth = 0.5) +
    ggplot2::geom_violin(width = 0.7, trim = FALSE, alpha = 0.35, linewidth = 0.9) +
    ggplot2::geom_hline(
      data = inter_method_line_data,
      ggplot2::aes(yintercept = yintercept, linetype = legend_label),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 1,
      show.legend = TRUE
    ) +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    ggplot2::geom_text(
      data = p_value_data,
      ggplot2::aes(x = method, y = y_pos, label = label),
      inherit.aes = FALSE,
      color = "black",
      vjust = 0,
      size = 3.6
    ) +
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
    ggplot2::scale_linetype_manual(
      values = stats::setNames("dashed", inter_method_legend_label),
      name = NULL
    ) +
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
      legend.position = "right",
      legend.justification = "center",
      legend.box.margin = ggplot2::margin(0, 0, 0, 8),
      plot.margin = ggplot2::margin(10, 24, 46, 10)
    )
}

silhouette_width_biplot <- function(meta1, meta2, name1, name2, for_pdf = FALSE) {
  required_meta_cols <- "silhouette_width"
  missing_meta1 <- setdiff(required_meta_cols, names(meta1))
  missing_meta2 <- setdiff(required_meta_cols, names(meta2))

  if (length(missing_meta1) > 0) {
    stop("Metadata for method `", name1, "` is missing required column `silhouette_width`.")
  }

  if (length(missing_meta2) > 0) {
    stop("Metadata for method `", name2, "` is missing required column `silhouette_width`.")
  }

  sil1 <- tibble::tibble(
    barcode = rownames(meta1),
    silhouette_1 = as.numeric(meta1$silhouette_width)
  )

  sil2 <- tibble::tibble(
    barcode = rownames(meta2),
    silhouette_2 = as.numeric(meta2$silhouette_width)
  )

  plot_df <- dplyr::full_join(sil1, sil2, by = "barcode") %>%
    dplyr::mutate(
      barcode_status = dplyr::case_when(
        !is.na(silhouette_1) & !is.na(silhouette_2) ~ "Shared barcode",
        !is.na(silhouette_1) & is.na(silhouette_2) ~ "Missing from one method",
        is.na(silhouette_1) & !is.na(silhouette_2) ~ "Missing from one method",
        TRUE ~ "Missing in both"
      ),
      silhouette_1 = dplyr::if_else(is.na(silhouette_1), -1, silhouette_1),
      silhouette_2 = dplyr::if_else(is.na(silhouette_2), -1, silhouette_2)
    ) %>%
    dplyr::filter(.data$barcode_status != "Missing in both")

  if (nrow(plot_df) == 0) {
    stop("No cell barcodes were available for the selected methods.")
  }

  status_levels <- c("Shared barcode", "Missing from one method")
  palette_values <- c(
    "Shared barcode" = "#1f78b4",
    "Missing from one method" = "grey60"
  )

  point_size <- if (for_pdf) 0.9 else 0.8
  alpha_value <- if (for_pdf) 0.22 else 0.18

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = silhouette_1,
      y = silhouette_2,
      color = factor(barcode_status, levels = status_levels)
    )
  ) +
    ggplot2::geom_vline(xintercept = -1, linetype = "dotted", color = "grey70", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -1, linetype = "dotted", color = "grey70", linewidth = 0.4) +
    ggplot2::geom_point(size = point_size, alpha = alpha_value) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey55", linewidth = 0.5) +
    ggplot2::coord_equal(xlim = c(-1, 1), ylim = c(-1, 1), expand = TRUE) +
    ggplot2::scale_color_manual(values = palette_values, drop = FALSE, name = NULL) +
    ggplot2::labs(
      title = "Cell-level silhouette width biplot",
      subtitle = paste0(scales::comma(nrow(plot_df)), " cell barcodes across two methods"),
      x = name1,
      y = name2
    ) +
    ggplot2::theme_bw(base_size = if (for_pdf) 12 else 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = ggplot2::element_text(size = if (for_pdf) 10 else 9),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(size = if (for_pdf) 2.8 else 2.4, alpha = 0.9)
      )
    )
}

build_silhouette_nj_tree <- function(meta_list, method_names, bootstrap = FALSE, n_boot = 100) {
  if (!is.list(meta_list) || length(meta_list) < 2) {
    stop("Provide at least two metadata tables to build a silhouette-width NJ tree.")
  }

  if (length(meta_list) != length(method_names)) {
    stop("`meta_list` and `method_names` must have the same length.")
  }

  if (length(unique(method_names)) != length(method_names)) {
    stop("`method_names` must be unique.")
  }

  if (!is.logical(bootstrap) || length(bootstrap) != 1 || is.na(bootstrap)) {
    stop("`bootstrap` must be TRUE or FALSE.")
  }

  if (!is.numeric(n_boot) || length(n_boot) != 1 || is.na(n_boot) || n_boot < 1) {
    stop("`n_boot` must be a single number greater than or equal to 1.")
  }

  n_boot <- as.integer(n_boot)

  missing_cols <- purrr::map_lgl(meta_list, ~ !("silhouette_width" %in% names(.x)))
  if (any(missing_cols)) {
    stop(
      "These methods are missing `silhouette_width` in metadata: ",
      paste(method_names[missing_cols], collapse = ", ")
    )
  }

  barcode_sets <- purrr::map(meta_list, rownames)
  shared_barcodes <- Reduce(intersect, barcode_sets)

  if (length(shared_barcodes) == 0) {
    stop("No shared barcodes were found across the selected methods.")
  }

  silhouette_mat <- vapply(meta_list, function(meta_tbl) {
    values <- as.numeric(meta_tbl[shared_barcodes, "silhouette_width", drop = TRUE])
    values
  }, numeric(length(shared_barcodes)))

  if (is.null(dim(silhouette_mat))) {
    silhouette_mat <- matrix(silhouette_mat, ncol = length(meta_list))
  }

  colnames(silhouette_mat) <- method_names
  rownames(silhouette_mat) <- shared_barcodes

  complete_rows <- stats::complete.cases(silhouette_mat)
  silhouette_mat <- silhouette_mat[complete_rows, , drop = FALSE]

  if (nrow(silhouette_mat) == 0) {
    stop("No shared barcodes had non-missing silhouette widths across all selected methods.")
  }

  rmsd <- function(x, y) {
    sqrt(mean((x - y)^2))
  }

  build_dist_from_matrix <- function(value_mat) {
    n_methods <- ncol(value_mat)
    dist_mat <- matrix(
      0,
      nrow = n_methods,
      ncol = n_methods,
      dimnames = list(method_names, method_names)
    )

    for (i in seq_len(n_methods)) {
      for (j in seq_len(i - 1)) {
        d <- rmsd(value_mat[, i], value_mat[, j])
        dist_mat[i, j] <- d
        dist_mat[j, i] <- d
      }
    }

    stats::as.dist(dist_mat)
  }

  dist_obj <- build_dist_from_matrix(silhouette_mat)
  nj_tree <- ape::nj(dist_obj)

  bootstrap_support <- NULL
  bootstrap_trees <- NULL
  if (bootstrap) {
    bootstrap_trees <- replicate(n_boot, {
      sampled_idx <- sample(seq_len(nrow(silhouette_mat)), size = nrow(silhouette_mat), replace = TRUE)
      sampled_mat <- silhouette_mat[sampled_idx, , drop = FALSE]
      ape::nj(build_dist_from_matrix(sampled_mat))
    }, simplify = FALSE)

    bootstrap_support <- ape::prop.clades(nj_tree, bootstrap_trees) / n_boot * 100
  }

  list(
    dist_matrix = dist_obj,
    nj_tree = nj_tree,
    bootstrap = bootstrap,
    n_boot = if (bootstrap) n_boot else 0L,
    bootstrap_support = bootstrap_support,
    bootstrap_trees = bootstrap_trees,
    shared_barcodes = shared_barcodes,
    shared_barcodes_used = rownames(silhouette_mat),
    dropped_shared_barcodes = setdiff(shared_barcodes, rownames(silhouette_mat))
  )
}

plot_silhouette_nj_tree <- function(nj_tree, bootstrap_support = NULL, for_pdf = FALSE) {
  ape::plot.phylo(
    nj_tree,
    main = "Neighbor-joining tree from silhouette width RMSD",
    cex = if (for_pdf) 0.78 else 0.72,
    font = 2,
    no.margin = FALSE
  )

  if (!is.null(bootstrap_support)) {
    ape::nodelabels(
      text = sprintf("%d", round(bootstrap_support)),
      frame = "none",
      adj = c(1.15, -0.2),
      cex = if (for_pdf) 0.85 else 0.75,
      col = "black"
    )
  }
}

ClusterStabilityVsSilhouettePlot <- function(meta_list, downsamp_list, method_labels,
                                             x_metric = "median_silhouette",
                                             y_metric = "median_max_jaccard",
                                             color_metric = "cluster_size",
                                             for_pdf = FALSE) {
  if (length(meta_list) != length(downsamp_list) || length(meta_list) != length(method_labels)) {
    stop("meta_list, downsamp_list, and method_labels must have the same length.")
  }

  if (length(meta_list) == 0) {
    stop("Provide at least one method.")
  }

  metric_specs <- list(
    median_silhouette = list(
      label = "Median silhouette width",
      limits = c(-1, 1),
      midpoint = 0
    ),
    median_max_jaccard = list(
      label = as.expression(expression(paste("Cluster ", stability[Jaccard]))),
      limits = c(0, 1),
      midpoint = 0.5
    ),
    cluster_size = list(
      label = "Cluster size",
      limits = NULL,
      midpoint = NULL
    )
  )

  selected_metrics <- c(x_metric, y_metric, color_metric)
  invalid_metrics <- setdiff(selected_metrics, names(metric_specs))
  if (length(invalid_metrics) > 0) {
    stop("Unknown silhouette/stability metric selection: ", paste(invalid_metrics, collapse = ", "))
  }

  if (anyDuplicated(selected_metrics) > 0) {
    stop("Choose three different metrics for x-axis, y-axis, and color.")
  }

  summary_df <- purrr::map_dfr(seq_along(meta_list), function(i) {
    meta_tbl <- meta_list[[i]]
    downsamp_tbl <- downsamp_list[[i]]
    method <- method_labels[[i]]
    required_meta_cols <- c("seurat_clusters", "silhouette_width")
    missing_meta_cols <- setdiff(required_meta_cols, names(meta_tbl))
    required_downsamp_cols <- c("clusterid", "max_jaccard")
    missing_downsamp_cols <- setdiff(required_downsamp_cols, names(downsamp_tbl))

    if (length(missing_meta_cols) > 0) {
      stop(
        "Metadata for method `",
        method,
        "` is missing required column",
        if (length(missing_meta_cols) > 1) "s: " else ": ",
        paste0("`", missing_meta_cols, "`", collapse = ", "),
        "."
      )
    }

    if (length(missing_downsamp_cols) > 0) {
      stop(
        "Downsampling summary for method `",
        method,
        "` is missing required column",
        if (length(missing_downsamp_cols) > 1) "s: " else ": ",
        paste0("`", missing_downsamp_cols, "`", collapse = ", "),
        "."
      )
    }

    size_summary <- meta_tbl %>%
      dplyr::group_by(clusterid = .data[["seurat_clusters"]]) %>%
      dplyr::summarise(
        cluster_size = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(clusterid = as.character(clusterid))

    sil_summary <- meta_tbl %>%
      dplyr::group_by(clusterid = .data[["seurat_clusters"]]) %>%
      dplyr::summarise(
        median_silhouette = median(.data[["silhouette_width"]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(clusterid = as.character(clusterid))

    jac_summary <- downsamp_tbl %>%
      dplyr::mutate(clusterid = as.character(clusterid)) %>%
      dplyr::group_by(clusterid) %>%
      dplyr::summarise(
        median_max_jaccard = median(max_jaccard, na.rm = TRUE),
        .groups = "drop"
      )

    sil_summary %>%
      dplyr::inner_join(jac_summary, by = "clusterid") %>%
      dplyr::inner_join(size_summary, by = "clusterid") %>%
      dplyr::mutate(method = method)
  })

  if (nrow(summary_df) == 0) {
    stop("No clusters remained after joining silhouette and Jaccard summaries.")
  }

  summary_df <- summary_df %>%
    dplyr::mutate(
      clusterid = factor(clusterid),
      method = factor(method, levels = method_labels)
    )

  x_spec <- metric_specs[[x_metric]]
  y_spec <- metric_specs[[y_metric]]
  color_spec <- metric_specs[[color_metric]]

  render_method_strip_label <- function(label, target_width, min_font_size, base_font_size, max_lines = 3) {
    label <- stringr::str_squish(as.character(label))

    if (!nzchar(label)) {
      return(label)
    }

    split_for_strip <- function(text) {
      text %>%
        stringr::str_replace_all("([+_-])", "\\1 ") %>%
        stringr::str_squish() %>%
        stringr::str_split("\\s+", simplify = FALSE) %>%
        purrr::pluck(1)
    }

    single_line_size <- base_font_size * target_width / max(1, nchar(label))
    tokens <- split_for_strip(label)
    has_break_opportunity <- length(tokens) > 1

    should_keep_single_line <- single_line_size >= min_font_size && nchar(label) <= target_width * 0.9
    if (should_keep_single_line) {
      return(label)
    }

    if (!has_break_opportunity) {
      return(label)
    }

    token_count <- length(tokens)
    max_breaks <- min(max_lines - 1, token_count - 1)

    if (max_breaks <= 0) {
      return(label)
    }

    best_label <- label
    best_score <- Inf

    for (break_count in seq_len(max_breaks)) {
      break_sets <- combn(token_count - 1, break_count, simplify = FALSE)

      for (breaks in break_sets) {
        starts <- c(1, breaks + 1)
        ends <- c(breaks, token_count)
        lines <- purrr::map2_chr(starts, ends, ~ paste(tokens[.x:.y], collapse = " ")) %>%
          stringr::str_replace_all("\\s+([+_-])$", "\\1")
        longest_line <- max(nchar(lines), na.rm = TRUE)
        wrapped_size <- base_font_size * target_width / max(1, longest_line)

        if (wrapped_size < min_font_size && wrapped_size <= single_line_size) {
          next
        }

        score <- sum((target_width - nchar(lines))^2) + break_count * 8 - wrapped_size * 6

        if (score < best_score) {
          best_score <- score
          best_label <- paste(lines, collapse = "\n")
        }
      }
    }

    best_label
  }

  mid_cs <- stats::median(summary_df$cluster_size, na.rm = TRUE)
  n_methods <- length(method_labels)
  facet_cols <- max(1, ceiling(sqrt(n_methods)))
  base_strip_size <- if (for_pdf) 11 else 10
  min_strip_size <- if (for_pdf) 8.5 else 8
  target_line_chars <- max(22, floor(if (for_pdf) 56 / facet_cols else 50 / facet_cols))

  label_line_width <- function(label) {
    lines <- stringr::str_split(as.character(label), "\n", simplify = FALSE)[[1]]
    max(nchar(lines), na.rm = TRUE)
  }

  candidate_labels <- purrr::map(
    1:3,
    ~ purrr::map_chr(
      method_labels,
      render_method_strip_label,
      target_width = target_line_chars,
      min_font_size = min_strip_size,
      base_font_size = base_strip_size,
      max_lines = .x
    )
  )
  candidate_longest_lines <- purrr::map_dbl(candidate_labels, ~ max(purrr::map_int(.x, label_line_width), na.rm = TRUE))
  candidate_sizes <- purrr::map_dbl(
    candidate_longest_lines,
    ~ min(base_strip_size, base_strip_size * target_line_chars / max(1, .x))
  )
  chosen_idx <- which(candidate_sizes >= min_strip_size)[1]
  allow_below_min_strip_size <- is.na(chosen_idx)
  if (is.na(chosen_idx)) {
    chosen_idx <- which.max(candidate_sizes)
  }

  rendered_method_labels <- candidate_labels[[chosen_idx]]
  longest_rendered_line <- candidate_longest_lines[[chosen_idx]]
  strip_text_size <- min(base_strip_size, base_strip_size * target_line_chars / max(1, longest_rendered_line))
  if (!allow_below_min_strip_size) {
    strip_text_size <- max(min_strip_size, strip_text_size)
  }

  summary_df <- summary_df %>%
    dplyr::mutate(
      method_display = factor(
        rendered_method_labels[match(as.character(method), method_labels)],
        levels = rendered_method_labels
      )
    )

  ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x = .data[[x_metric]],
      y = .data[[y_metric]],
      fill = .data[[color_metric]]
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      color = "black",
      size = if (for_pdf) 2.6 else 2.2,
      alpha = 1
    ) +
    ggplot2::facet_wrap(
      ~ method_display,
      scales = "fixed"
    ) +
    ggplot2::scale_x_continuous(limits = x_spec$limits) +
    ggplot2::scale_y_continuous(limits = y_spec$limits) +
    ggplot2::scale_fill_gradient2(
      name = color_spec$label,
      low = "blue",
      mid = "white",
      high = "firebrick",
      midpoint = if (is.null(color_spec$midpoint)) stats::median(summary_df[[color_metric]], na.rm = TRUE) else color_spec$midpoint
    ) +
    ggplot2::xlab(x_spec$label) +
    ggplot2::ylab(y_spec$label) +
    ggplot2::theme_bw(base_size = if (for_pdf) 12 else 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = strip_text_size,
        hjust = 0.5,
        vjust = 0.5,
        lineheight = 0.95,
        margin = ggplot2::margin(5, 0, 5, 0)
      ),
      legend.position = "right"
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
  max_method_label_chars <- max(nchar(method_levels), 1)

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

  y_axis_text_size <- if (for_pdf) {
    max(6.2, 10 - 0.13 * max(max_method_label_chars - 22, 0))
  } else {
    max(5.5, 8 - 0.11 * max(max_method_label_chars - 22, 0))
  }
  plot_title_size <- if (for_pdf) 13 else 11
  tile_height <- if (for_pdf) 0.84 else 0.92
  tile_width <- 1
  top_margin <- if (for_pdf) 28 else 48
  bottom_margin <- if (identical(sort_mode, "cluster") && nrow(cluster_annotations) > 0) {
    if (for_pdf) 56 else 84
  } else {
    if (for_pdf) 26 else 42
  }
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
      axis.text.y = ggplot2::element_text(size = y_axis_text_size, color = "black", margin = ggplot2::margin(r = 4)),
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
