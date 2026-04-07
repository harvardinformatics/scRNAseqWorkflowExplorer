library(shiny)
library(tidyverse)
library(Seurat)
library(ggrepel)

jaccard_similarity <- function(set1, set2) {
  intersect_length <- length(intersect(set1, set2))
  union_length <- length(set1) + length(set2) - intersect_length
  intersect_length / union_length
}

standardize_bootstrap_cols <- function(df) {
  nm <- names(df)
  nm_std <- nm %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")

  names(df) <- nm_std

  if ("maxjaccard" %in% names(df) && !("max_jaccard" %in% names(df))) {
    df <- dplyr::rename(df, max_jaccard = maxjaccard)
  }

  if ("boostrap_number" %in% names(df) && !("bootstrap_number" %in% names(df))) {
    df <- dplyr::rename(df, bootstrap_number = boostrap_number)
  }

  if ("boostrapnumber" %in% names(df) && !("bootstrap_number" %in% names(df))) {
    df <- dplyr::rename(df, bootstrap_number = boostrapnumber)
  }

  if ("bootstrapnumber" %in% names(df) && !("bootstrap_number" %in% names(df))) {
    df <- dplyr::rename(df, bootstrap_number = bootstrapnumber)
  }

  required <- c("clusterid", "max_jaccard", "bootstrap_number")
  missing <- setdiff(required, names(df))

  if (length(missing) > 0) {
    stop(
      "Bootstrap TSV is missing required columns after normalization: ",
      paste(missing, collapse = ", ")
    )
  }

  df
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

find_method_pairs <- function(data_dir = "data") {
  rds_paths <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)

  if (length(rds_paths) == 0) {
    return(tibble::tibble(label = character(), rds = character(), boot = character()))
  }

  tibble::tibble(rds = rds_paths) %>%
    dplyr::mutate(
      boot = purrr::map_chr(rds, infer_bootstrap_path),
      exists_boot = file.exists(boot),
      label = basename(rds)
    ) %>%
    dplyr::filter(exists_boot) %>%
    dplyr::select(label, rds, boot) %>%
    dplyr::arrange(label)
}

infer_bootstrap_path <- function(rds_path) {
  base_no_ext <- stringr::str_replace(rds_path, "\\.rds$", "")
  candidates <- c(
    paste0(base_no_ext, "_bootstraps.tsv"),
    paste0(base_no_ext, "_clusterbootstraps.tsv")
  )

  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    return(candidates[[1]])
  }
  existing[[1]]
}

read_seurat_meta <- function(rds_path) {
  obj <- readRDS(rds_path)

  if (!inherits(obj, "Seurat")) {
    stop("RDS is not a Seurat object: ", basename(rds_path))
  }

  if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
    stop("Missing `seurat_clusters` in Seurat metadata: ", basename(rds_path))
  }

  obj@meta.data
}

read_bootstrap_tsv <- function(tsv_path) {
  readr::read_tsv(tsv_path, show_col_types = FALSE) %>%
    standardize_bootstrap_cols()
}

extract_feature_table <- function(obj, assay_name) {
  assay_obj <- obj[[assay_name]]
  feature_ids <- rownames(SeuratObject::GetAssayData(obj, assay = assay_name, layer = "counts"))
  assay_symbols <- row.names(assay_obj)

  gene_symbols <- if (
    !is.null(assay_symbols) &&
    length(assay_symbols) == length(feature_ids)
  ) {
    assay_symbols
  } else {
    feature_ids
  }

  tibble::tibble(
    feature_id = feature_ids,
    gene_symbol = gene_symbols
  )
}

select_expression_assay <- function(obj) {
  assay_names <- names(obj@assays)

  if ("RNA" %in% assay_names) {
    return("RNA")
  }

  SeuratObject::DefaultAssay(obj)
}

read_method_bundle <- function(rds_path, boot_path) {
  obj <- readRDS(rds_path)

  if (!inherits(obj, "Seurat")) {
    stop("RDS is not a Seurat object: ", basename(rds_path))
  }

  if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
    stop("Missing `seurat_clusters` in Seurat metadata: ", basename(rds_path))
  }

  assay_name <- select_expression_assay(obj)
  counts <- SeuratObject::GetAssayData(obj, assay = assay_name, layer = "counts")
  lib_size <- Matrix::colSums(counts)
  feature_table <- extract_feature_table(obj, assay_name)

  list(
    label = stringr::str_remove(basename(rds_path), "\\.rds$"),
    rds = rds_path,
    boot = boot_path,
    assay = assay_name,
    meta = obj@meta.data,
    bootstraps = read_bootstrap_tsv(boot_path),
    counts = counts,
    lib_size = lib_size,
    barcodes = colnames(counts),
    features = rownames(counts),
    feature_table = feature_table
  )
}

resolve_feature_id <- function(method, gene_symbol) {
  match_tbl <- method$feature_table %>%
    dplyr::filter(.data$gene_symbol == .env$gene_symbol)

  if (nrow(match_tbl) == 0) {
    return(NULL)
  }

  match_tbl$feature_id[[1]]
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

move_sort_method_to_bottom <- function(method_levels, sort_method) {
  if (is.null(sort_method) || !nzchar(sort_method) || !(sort_method %in% method_levels)) {
    return(method_levels)
  }

  c(setdiff(method_levels, sort_method), sort_method)
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

make_expression_heatmap_data <- function(method_data, gene_symbol, sort_method = NULL, sort_mode = "expression", min_cluster_size = 1) {
  barcode_union <- sort(unique(unlist(purrr::map(method_data, "barcodes"))))
  method_levels <- move_sort_method_to_bottom(purrr::map_chr(method_data, "label"), sort_method)
  cluster_gap <- 0

  plot_data <- purrr::map_dfr(method_data, function(method) {
    values <- rep(NA_real_, length(barcode_union))
    names(values) <- barcode_union
    barcode_present <- barcode_union %in% method$barcodes
    cell_status <- ifelse(barcode_present, "expression_missing", "barcode_missing")

    feature_id <- resolve_feature_id(method, gene_symbol)

    if (!is.null(feature_id)) {
      gene_counts <- as.numeric(method$counts[feature_id, , drop = TRUE])
      normalized <- log1p((gene_counts / pmax(method$lib_size, 1)) * 10000)
      values[method$barcodes] <- normalized
      cell_status[barcode_present] <- "expression_present"
    }

    tibble::tibble(
      method = method$label,
      barcode = factor(barcode_union, levels = barcode_union),
      expression = values,
      cell_status = factor(
        cell_status,
        levels = c("barcode_missing", "expression_missing", "expression_present")
      )
    )
  }) %>%
    dplyr::mutate(
      method = factor(method, levels = method_levels)
    )

  cluster_annotations <- tibble::tibble()
  cluster_boundaries <- numeric()
  excluded_clusters <- tibble::tibble()

  if (!is.null(sort_method) && nzchar(sort_method) && sort_method %in% levels(plot_data$method)) {
    sort_method_data <- method_data[[match(sort_method, purrr::map_chr(method_data, "label"))]]
    sort_reference <- plot_data %>%
      dplyr::filter(method == sort_method) %>%
      dplyr::mutate(barcode = as.character(barcode))

    if (identical(sort_mode, "cluster")) {
      cluster_lookup <- tibble::tibble(
        barcode = rownames(sort_method_data$meta),
        cluster = as.character(sort_method_data$meta$seurat_clusters)
      )

      cluster_sizes <- cluster_lookup %>%
        dplyr::filter(!is.na(cluster)) %>%
        dplyr::count(cluster, name = "cluster_size") %>%
        dplyr::arrange(dplyr::desc(cluster_size), cluster)

      included_clusters <- cluster_sizes %>%
        dplyr::filter(cluster_size >= min_cluster_size) %>%
        dplyr::pull(cluster)

      excluded_clusters <- cluster_sizes %>%
        dplyr::filter(cluster_size < min_cluster_size) %>%
        dplyr::arrange(cluster_size, cluster)

      sort_reference <- sort_reference %>%
        dplyr::left_join(cluster_lookup, by = "barcode") %>%
        dplyr::mutate(
          cluster_original = cluster,
          cluster = dplyr::if_else(cluster %in% included_clusters, cluster, NA_character_),
          cluster = factor(cluster, levels = arrange_cluster_levels(included_clusters))
        )

      cluster_summary <- sort_reference %>%
        dplyr::filter(cell_status == "expression_present", !is.na(cluster)) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(cluster_median = stats::median(expression, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(cluster_median, cluster)

      cluster_levels_present <- cluster_summary %>%
        dplyr::pull(cluster) %>%
        as.character()

      sort_reference <- sort_reference %>%
        dplyr::left_join(cluster_summary, by = "cluster") %>%
        dplyr::mutate(
          missing_rank = dplyr::case_when(
            cell_status == "barcode_missing" ~ 2L,
            cell_status == "expression_missing" ~ 1L,
            TRUE ~ 0L
          ),
          expression_within_cluster = dplyr::if_else(
            missing_rank == 0L,
            expression,
            Inf
          )
        )

      clustered_reference <- sort_reference %>%
        dplyr::filter(!is.na(cluster)) %>%
        dplyr::mutate(
          cluster = factor(as.character(cluster), levels = cluster_levels_present)
        ) %>%
        dplyr::arrange(cluster, missing_rank, expression_within_cluster, barcode) %>%
        dplyr::mutate(
          cluster_break = dplyr::if_else(
            dplyr::row_number() == 1L,
            0,
            dplyr::if_else(as.character(cluster) != dplyr::lag(as.character(cluster)), 1, 0)
          ),
          gap_index = cumsum(cluster_break),
          display_x = dplyr::row_number() + gap_index * cluster_gap
        )

      unclustered_reference <- sort_reference %>%
        dplyr::filter(is.na(cluster_original)) %>%
        dplyr::arrange(missing_rank, expression, barcode)

      if (nrow(unclustered_reference) > 0) {
        unclustered_start <- if (nrow(clustered_reference) > 0) {
          max(clustered_reference$display_x)
        } else {
          0
        }

        unclustered_reference <- unclustered_reference %>%
          dplyr::mutate(display_x = unclustered_start + dplyr::row_number())
      }

      sort_reference <- dplyr::bind_rows(clustered_reference, unclustered_reference)

      cluster_annotations <- clustered_reference %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(
          x = mean(range(display_x)),
          x_start = min(display_x),
          x_end = max(display_x),
          xmin = min(display_x) - 0.5,
          xmax = max(display_x) + 0.5,
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          label = as.character(cluster),
          sort_method = sort_method
        )

      cluster_boundaries <- clustered_reference %>%
        dplyr::arrange(display_x) %>%
        dplyr::mutate(
          previous_cluster = dplyr::lag(as.character(cluster)),
          previous_x = dplyr::lag(display_x),
          boundary_x = (previous_x + display_x) / 2
        ) %>%
        dplyr::filter(!is.na(previous_cluster), as.character(cluster) != previous_cluster) %>%
        dplyr::transmute(boundary_x = floor(boundary_x) + 0.5) %>%
        dplyr::distinct(boundary_x) %>%
        dplyr::pull(boundary_x)

      if (nrow(cluster_annotations) > 0) {
        cluster_palette <- make_contrast_cluster_palette(nrow(cluster_annotations))
        cluster_annotations <- cluster_annotations %>%
          dplyr::arrange(x) %>%
          dplyr::mutate(
            cluster_color = cluster_palette,
            label_y = dplyr::if_else(dplyr::row_number() %% 2 == 1, -0.24, -0.48)
          )
      }
    } else {
      sort_reference <- sort_reference %>%
        dplyr::mutate(
          missing_rank = ifelse(is.na(expression), 1L, 0L)
        ) %>%
        dplyr::arrange(missing_rank, expression, barcode) %>%
        dplyr::mutate(display_x = dplyr::row_number())
    }

    barcode_levels <- sort_reference %>%
      dplyr::pull(barcode) %>%
      as.character()

    barcode_positions <- sort_reference %>%
      dplyr::transmute(
        barcode = as.character(barcode),
        display_x = display_x
      )

    plot_data <- plot_data %>%
      dplyr::mutate(
        barcode = factor(as.character(barcode), levels = barcode_levels),
        display_x = barcode_positions$display_x[match(as.character(barcode), barcode_positions$barcode)]
      ) %>%
      dplyr::filter(!is.na(display_x))
  }

  list(
    plot_data = plot_data,
    cluster_annotations = cluster_annotations,
    cluster_boundaries = cluster_boundaries,
    excluded_clusters = excluded_clusters
  )
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
  # Keep tiles shorter than the row spacing so each method is visually separated.
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

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .selectize-control.single .selectize-input,
      .selectize-dropdown .option {
        font-size: clamp(9px, 0.95vw, 12px);
      }
      .selectize-control.single .selectize-input > div.item {
        font-size: clamp(9px, 0.95vw, 12px);
      }
    "))
  ),
  titlePanel("scRNA-seq method comparison"),
  tabsetPanel(
    id = "analysis_tabs",
    tabPanel(
      "Cluster stability",
      sidebarLayout(
        sidebarPanel(
          width = 5,
          helpText("Select two methods (.rds) with matching bootstrap TSV files."),
          selectInput("method1", "Method 1", choices = NULL),
          selectInput("method2", "Method 2", choices = NULL),
          numericInput("threshold", "Minimum Jaccard threshold", value = 0.6, min = 0, max = 1, step = 0.01),
          actionButton("refresh", "Refresh file list")
        ),
        mainPanel(
          width = 7,
          plotOutput("stability_plot", height = "500px"),
          verbatimTextOutput("status")
        )
      )
    ),
    tabPanel(
      "Gene heatmap",
      sidebarLayout(
        sidebarPanel(
          width = 4,
          helpText("Choose a gene to compare log-normalized expression across all methods and the union of all cell barcodes."),
          shiny::tagAppendAttributes(
            textInput("heatmap_gene", "Gene symbol", value = ""),
            list = "gene-symbol-options",
            autocomplete = "off"
          ),
          uiOutput("gene_symbol_datalist"),
          selectInput(
            "heatmap_sort_mode",
            "Sort mode",
            choices = c(
              "Expression" = "expression",
              "Cluster median expression" = "cluster"
            ),
            selected = "expression"
          ),
          selectInput("heatmap_sort_method", "Sort barcodes by method", choices = c()),
          conditionalPanel(
            condition = "input.heatmap_sort_mode === 'cluster'",
            numericInput("heatmap_min_cluster_size", "Minimum cells per cluster", value = 1, min = 1, step = 1)
          ),
          actionButton("refresh_heatmap", "Refresh gene list")
          ,
          downloadButton("download_heatmap_pdf", "Download hi-res PDF")
        ),
        mainPanel(
          width = 8,
          tags$div(
            style = "position: relative;",
            plotOutput(
              "expression_heatmap",
              height = "590px",
              hover = hoverOpts("expression_heatmap_hover", delay = 80, delayType = "debounce")
            ),
            uiOutput("cluster_hover_tooltip")
          ),
          conditionalPanel(
            condition = "input.heatmap_sort_mode === 'cluster'",
            tags$div(
              style = "margin-top: 4px; margin-bottom: 10px; font-size: 12px; color: #555;",
              "Hover over cluster annotations to obtain the cluster id"
            )
          ),
          uiOutput("excluded_clusters_header"),
          tableOutput("excluded_clusters_table"),
          verbatimTextOutput("heatmap_status")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  method_pairs <- reactiveVal(find_method_pairs("data"))

  refresh_choices <- function() {
    pairs <- find_method_pairs("data")
    method_pairs(pairs)

    if (nrow(pairs) > 0) {
      choice_map <- stats::setNames(pairs$rds, pairs$label)
      current_method1 <- isolate(input$method1)
      current_method2 <- isolate(input$method2)

      selected_method1 <- if (!is.null(current_method1) && current_method1 %in% pairs$rds) {
        current_method1
      } else {
        pairs$rds[[1]]
      }

      fallback_method2 <- if (nrow(pairs) >= 2) pairs$rds[[2]] else pairs$rds[[1]]
      selected_method2 <- if (!is.null(current_method2) && current_method2 %in% pairs$rds) {
        current_method2
      } else {
        fallback_method2
      }

      if (identical(selected_method1, selected_method2) && nrow(pairs) >= 2) {
        selected_method2 <- pairs$rds[[which(pairs$rds != selected_method1)[1]]]
      }

      updateSelectInput(session, "method1", choices = choice_map, selected = selected_method1)
      updateSelectInput(session, "method2", choices = choice_map, selected = selected_method2)
    } else {
      updateSelectInput(session, "method1", choices = c())
      updateSelectInput(session, "method2", choices = c())
    }
  }

  observeEvent(input$refresh, refresh_choices(), ignoreInit = TRUE)
  observeEvent(input$refresh_heatmap, {
    current_gene <- isolate(input$heatmap_gene)
    method_pairs(find_method_pairs("data"))

    if (!is.null(current_gene) && nzchar(current_gene)) {
      updateTextInput(session, "heatmap_gene", value = current_gene)
    }
  }, ignoreInit = TRUE)
  observe(refresh_choices())

  loaded_data <- reactive({
    req(input$method1, input$method2)

    validate(
      need(input$method1 != input$method2, "Pick two different methods."),
      need(file.exists(input$method1), "Method 1 file does not exist."),
      need(file.exists(input$method2), "Method 2 file does not exist.")
    )

    boot1 <- infer_bootstrap_path(input$method1)
    boot2 <- infer_bootstrap_path(input$method2)

    validate(
      need(file.exists(boot1), paste("Missing bootstrap TSV for method 1:", basename(boot1))),
      need(file.exists(boot2), paste("Missing bootstrap TSV for method 2:", basename(boot2)))
    )

    list(
      meta1 = read_seurat_meta(input$method1),
      meta2 = read_seurat_meta(input$method2),
      b1 = read_bootstrap_tsv(boot1),
      b2 = read_bootstrap_tsv(boot2),
      l1 = stringr::str_remove(basename(input$method1), "\\.rds$"),
      l2 = stringr::str_remove(basename(input$method2), "\\.rds$")
    )
  })

  all_method_data <- reactive({
    pairs <- method_pairs()

    validate(
      need(nrow(pairs) > 0, "No valid method pairs found in ./data.")
    )

    purrr::map2(pairs$rds, pairs$boot, read_method_bundle)
  })

  output$gene_symbol_datalist <- renderUI({
    methods <- all_method_data()
    gene_choices <- methods %>%
      purrr::map("feature_table") %>%
      dplyr::bind_rows() %>%
      dplyr::distinct(gene_symbol) %>%
      dplyr::arrange(gene_symbol) %>%
      dplyr::pull(gene_symbol)

    tags$datalist(
      id = "gene-symbol-options",
      lapply(gene_choices, function(gene) {
        tags$option(value = gene)
      })
    )
  })

  observe({
    methods <- all_method_data()
    method_labels <- purrr::map_chr(methods, "label")
    current_sort <- isolate(input$heatmap_sort_method)
    selected_sort <- if (!is.null(current_sort) && current_sort %in% method_labels) {
      current_sort
    } else {
      method_labels[[1]]
    }

    updateSelectInput(
      session,
      "heatmap_sort_method",
      choices = stats::setNames(method_labels, method_labels),
      selected = selected_sort
    )
  })

  heatmap_data <- reactive({
    methods <- all_method_data()
    req(input$heatmap_gene)
    validate(
      need(nzchar(input$heatmap_gene), "Type a gene symbol to draw the heatmap.")
    )

    make_expression_heatmap_data(
      methods,
      input$heatmap_gene,
      input$heatmap_sort_method,
      input$heatmap_sort_mode,
      min_cluster_size = if (is.null(input$heatmap_min_cluster_size)) 1 else input$heatmap_min_cluster_size
    )
  })

  output$stability_plot <- renderPlot({
    dat <- loaded_data()

    MakeInterVsIntraStablePlot(
      meta1 = dat$meta1,
      meta2 = dat$meta2,
      bootstraps1 = dat$b1,
      bootstraps2 = dat$b2,
      threshold = input$threshold,
      label1 = dat$l1,
      label2 = dat$l2
    )
  })

  output$expression_heatmap <- renderPlot({
    heatmap_obj <- heatmap_data()
    plot_data <- heatmap_obj$plot_data
    validate(
      need(nrow(plot_data) > 0, "No heatmap data available for the selected gene."),
      need(any(!is.na(plot_data$expression)), "Selected gene symbol was not found in the loaded methods.")
    )

    make_expression_heatmap_plot(
      plot_data,
      input$heatmap_gene,
      cluster_annotations = heatmap_obj$cluster_annotations,
      cluster_boundaries = heatmap_obj$cluster_boundaries,
      sort_mode = input$heatmap_sort_mode,
      sort_method = input$heatmap_sort_method
    )
  }, res = 110)

  hovered_cluster_info <- reactive({
    heatmap_obj <- heatmap_data()
    cluster_annotations <- heatmap_obj$cluster_annotations
    hover <- input$expression_heatmap_hover

    if (
      is.null(hover) ||
      !identical(input$heatmap_sort_mode, "cluster") ||
      nrow(cluster_annotations) == 0 ||
      is.null(hover$x) || is.null(hover$y)
    ) {
      return(NULL)
    }

    hovered_cluster <- cluster_annotations %>%
      dplyr::filter(
        .data$xmin <= hover$x,
        .data$xmax >= hover$x,
        hover$y >= -0.08,
        hover$y <= 0.36
      ) %>%
      dplyr::slice(1)

    if (nrow(hovered_cluster) == 0 || is.null(hover$coords_css$x) || is.null(hover$coords_css$y)) {
      return(NULL)
    }

    list(
      cluster = hovered_cluster,
      hover = hover
    )
  })

  output$cluster_hover_tooltip <- renderUI({
    hover_info <- hovered_cluster_info()

    if (is.null(hover_info)) {
      return(NULL)
    }

    tags$div(
      style = paste(
        "position:absolute;",
        sprintf("left:%spx;", hover_info$hover$coords_css$x + 12),
        sprintf("top:%spx;", hover_info$hover$coords_css$y - 12),
        "pointer-events:none;",
        "background:rgba(255,255,255,0.96);",
        "border:1px solid #bdbdbd;",
        "border-radius:4px;",
        "padding:6px 8px 8px 8px;",
        "font-size:12px;",
        "width:320px;",
        "box-shadow:0 2px 8px rgba(0,0,0,0.12);",
        sep = " "
      ),
      plotOutput("cluster_hover_preview", width = "300px", height = "170px")
    )
  })

  output$cluster_hover_preview <- renderPlot({
    hover_info <- hovered_cluster_info()
    req(hover_info)

    heatmap_obj <- heatmap_data()
    make_cluster_preview_plot(
      heatmap_obj$plot_data,
      hover_info$cluster,
      input$heatmap_gene,
      sort_method = input$heatmap_sort_method
    )
  }, res = 120)

  output$excluded_clusters_header <- renderUI({
    heatmap_obj <- heatmap_data()
    excluded_clusters <- heatmap_obj$excluded_clusters

    if (!identical(input$heatmap_sort_mode, "cluster") || nrow(excluded_clusters) == 0) {
      return(NULL)
    }

    tags$div(
      style = "margin-top: 8px; font-size: 13px;",
      tags$strong("Excluded clusters")
    )
  })

  output$excluded_clusters_table <- renderTable({
    heatmap_obj <- heatmap_data()
    excluded_clusters <- heatmap_obj$excluded_clusters

    if (!identical(input$heatmap_sort_mode, "cluster") || nrow(excluded_clusters) == 0) {
      return(NULL)
    }

    excluded_clusters %>%
      dplyr::rename(
        Cluster = cluster,
        `Cell barcodes` = cluster_size
      )
  }, striped = TRUE, bordered = TRUE, spacing = "xs")

  output$download_heatmap_pdf <- downloadHandler(
    filename = function() {
      gene_label <- if (!is.null(input$heatmap_gene) && nzchar(input$heatmap_gene)) {
        input$heatmap_gene
      } else {
        "heatmap"
      }
      paste0("gene-heatmap-", gene_label, ".pdf")
    },
    content = function(file) {
      heatmap_obj <- heatmap_data()
      plot_data <- heatmap_obj$plot_data
      validate(
        need(nrow(plot_data) > 0, "No heatmap data available for the selected gene."),
        need(any(!is.na(plot_data$expression)), "Selected gene symbol was not found in the loaded methods.")
      )

      n_methods <- length(levels(plot_data$method))
      pdf_width <- 10
      pdf_height <- max(6.5, min(9, n_methods * 0.55 + 3.5))

      grDevices::pdf(file, width = pdf_width, height = pdf_height, onefile = TRUE)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(make_expression_heatmap_plot(
        plot_data,
        input$heatmap_gene,
        cluster_annotations = heatmap_obj$cluster_annotations,
        cluster_boundaries = heatmap_obj$cluster_boundaries,
        sort_mode = input$heatmap_sort_mode,
        sort_method = input$heatmap_sort_method,
        for_pdf = TRUE
      ))
    }
  )

  output$status <- renderText({
    pairs <- method_pairs()

    if (nrow(pairs) == 0) {
      return("No valid method pairs found in ./data. Expected .rds plus matching bootstrap TSV file.")
    }

    paste0("Detected ", nrow(pairs), " method file pairs in ./data")
  })

  output$heatmap_status <- renderText({
    heatmap_obj <- heatmap_data()
    methods <- all_method_data()
    barcode_union <- sort(unique(unlist(purrr::map(methods, "barcodes"))))
    feature_union <- methods %>%
      purrr::map("feature_table") %>%
      dplyr::bind_rows() %>%
      dplyr::distinct(gene_symbol) %>%
      nrow()
    cluster_filter_text <- if (identical(input$heatmap_sort_mode, "cluster")) {
      paste0(
        " Cluster filter: minimum ",
        if (is.null(input$heatmap_min_cluster_size)) 1 else input$heatmap_min_cluster_size,
        " cells; excluded ",
        nrow(heatmap_obj$excluded_clusters),
        " clusters."
      )
    } else {
      ""
    }

    paste0(
      "Loaded ", length(methods), " methods. Heatmap rows are methods, columns are the union of ",
      length(barcode_union), " cell barcodes, and gene choices cover ",
      feature_union, " gene symbols.",
      cluster_filter_text
    )
  })
}

shinyApp(ui = ui, server = server)
