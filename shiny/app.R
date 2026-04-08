library(shiny)
library(tidyverse)
library(Seurat)
library(ggrepel)
source("functions.R")

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

move_sort_method_to_bottom <- function(method_levels, sort_method) {
  if (is.null(sort_method) || !nzchar(sort_method) || !(sort_method %in% method_levels)) {
    return(method_levels)
  }

  c(setdiff(method_levels, sort_method), sort_method)
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
        droplevels() %>%
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
          shiny::tagList(
            helpText("Select two Seurat objects representing different workflows (.rds) with matching clustering-from-downsampling (.tsv) files."),
            helpText("Based upon the user-specified Jaccard similarity threshold:"),
            tags$ul(
              tags$li("Distributions of the number of intra-method stable clusters are calculated, where stability is defined by similarity to clusters in the original data, and plotted as violin-boxplot combinations."),
              tags$li("From the original full-data objects, inter-method stable clusters are calculated in a similar fashion and are represented as a horizontal dashed line."),
              tags$li(
                HTML("For each method, the null hypothesis <i>H</i><sub>0</sub>: number of intra-method stable clusters = number of inter-method stable clusters is evaluated with a Wilcoxon signed-rank test.")
              )
            )
          ),
          selectInput("method1", "Method 1", choices = NULL),
          textInput("method1_label", "Method 1 label", value = "", placeholder = "Defaults to selected file name"),
          selectInput("method2", "Method 2", choices = NULL),
          textInput("method2_label", "Method 2 label", value = "", placeholder = "Defaults to selected file name"),
          numericInput("threshold", "Minimum Jaccard threshold", value = 0.6, min = 0, max = 1, step = 0.01),
          actionButton("refresh", "Refresh file list"),
          downloadButton("download_stability_pdf", "Download hi-res PDF")
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
    ),
    tabPanel(
      "Cluster similarity heatmap",
      sidebarLayout(
        sidebarPanel(
          width = 4,
          helpText("Select two methods and plot cluster-level Jaccard similarity between their Seurat clusterings."),
          selectInput("jaccard_method1", "Method 1", choices = NULL),
          textInput("jaccard_label1", "Method 1 axis label", value = ""),
          selectInput("jaccard_method2", "Method 2", choices = NULL),
          textInput("jaccard_label2", "Method 2 axis label", value = ""),
          numericInput(
            "jaccard_threshold",
            "Minimum Jaccard threshold for cell labels",
            value = 0.6,
            min = 0,
            max = 1,
            step = 0.01
          ),
          actionButton("refresh_jaccard", "Refresh method list"),
          downloadButton("download_jaccard_pdf", "Download hi-res PDF")
        ),
        mainPanel(
          width = 8,
          plotOutput("jaccard_heatmap", height = "620px"),
          verbatimTextOutput("jaccard_status")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  method_pairs <- reactiveVal(find_method_pairs("data"))
  jaccard_label_defaults <- reactiveValues(label1 = NULL, label2 = NULL)

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

      current_jaccard_method1 <- isolate(input$jaccard_method1)
      current_jaccard_method2 <- isolate(input$jaccard_method2)

      selected_jaccard_method1 <- if (!is.null(current_jaccard_method1) && current_jaccard_method1 %in% pairs$rds) {
        current_jaccard_method1
      } else {
        selected_method1
      }

      selected_jaccard_method2 <- if (!is.null(current_jaccard_method2) && current_jaccard_method2 %in% pairs$rds) {
        current_jaccard_method2
      } else {
        selected_method2
      }

      if (identical(selected_jaccard_method1, selected_jaccard_method2) && nrow(pairs) >= 2) {
        selected_jaccard_method2 <- pairs$rds[[which(pairs$rds != selected_jaccard_method1)[1]]]
      }

      updateSelectInput(session, "method1", choices = choice_map, selected = selected_method1)
      updateSelectInput(session, "method2", choices = choice_map, selected = selected_method2)
      updateSelectInput(session, "jaccard_method1", choices = choice_map, selected = selected_jaccard_method1)
      updateSelectInput(session, "jaccard_method2", choices = choice_map, selected = selected_jaccard_method2)
    } else {
      updateSelectInput(session, "method1", choices = c())
      updateSelectInput(session, "method2", choices = c())
      updateSelectInput(session, "jaccard_method1", choices = c())
      updateSelectInput(session, "jaccard_method2", choices = c())
    }
  }

  observeEvent(input$refresh, refresh_choices(), ignoreInit = TRUE)
  observeEvent(input$refresh_jaccard, refresh_choices(), ignoreInit = TRUE)
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
      b2 = read_bootstrap_tsv(boot2)
    )
  })

  stability_plot_labels <- reactive({
    req(input$method1, input$method2)

    list(
      l1 = if (!is.null(input$method1_label) && nzchar(trimws(input$method1_label))) {
        trimws(input$method1_label)
      } else {
        stringr::str_remove(basename(input$method1), "\\.rds$")
      },
      l2 = if (!is.null(input$method2_label) && nzchar(trimws(input$method2_label))) {
        trimws(input$method2_label)
      } else {
        stringr::str_remove(basename(input$method2), "\\.rds$")
      }
    )
  })

  loaded_jaccard_data <- reactive({
    req(input$jaccard_method1, input$jaccard_method2)

    validate(
      need(input$jaccard_method1 != input$jaccard_method2, "Pick two different methods."),
      need(file.exists(input$jaccard_method1), "Method 1 file does not exist."),
      need(file.exists(input$jaccard_method2), "Method 2 file does not exist.")
    )

    list(
      meta1 = read_seurat_meta(input$jaccard_method1),
      meta2 = read_seurat_meta(input$jaccard_method2),
      l1 = stringr::str_remove(basename(input$jaccard_method1), "\\.rds$"),
      l2 = stringr::str_remove(basename(input$jaccard_method2), "\\.rds$")
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

  observe({
    dat <- loaded_jaccard_data()

    current_label1 <- isolate(input$jaccard_label1)
    current_label2 <- isolate(input$jaccard_label2)

    if (
      is.null(current_label1) ||
      !nzchar(current_label1) ||
      identical(current_label1, jaccard_label_defaults$label1)
    ) {
      updateTextInput(session, "jaccard_label1", value = dat$l1)
    }

    if (
      is.null(current_label2) ||
      !nzchar(current_label2) ||
      identical(current_label2, jaccard_label_defaults$label2)
    ) {
      updateTextInput(session, "jaccard_label2", value = dat$l2)
    }

    jaccard_label_defaults$label1 <- dat$l1
    jaccard_label_defaults$label2 <- dat$l2
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
    labels <- stability_plot_labels()

    MakeInterVsIntraStablePlot(
      meta1 = dat$meta1,
      meta2 = dat$meta2,
      bootstraps1 = dat$b1,
      bootstraps2 = dat$b2,
      threshold = input$threshold,
      label1 = labels$l1,
      label2 = labels$l2
    )
  })

  output$download_stability_pdf <- downloadHandler(
    filename = function() {
      labels <- stability_plot_labels()
      label1_slug <- gsub("[^A-Za-z0-9_-]+", "-", labels$l1)
      label2_slug <- gsub("[^A-Za-z0-9_-]+", "-", labels$l2)
      paste0(
        "cluster-stability-",
        label1_slug,
        "-vs-",
        label2_slug,
        ".pdf"
      )
    },
    content = function(file) {
      dat <- loaded_data()
      labels <- stability_plot_labels()
      grDevices::pdf(file, width = 9.5, height = 6.5, onefile = TRUE)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(
        MakeInterVsIntraStablePlot(
          meta1 = dat$meta1,
          meta2 = dat$meta2,
          bootstraps1 = dat$b1,
          bootstraps2 = dat$b2,
          threshold = input$threshold,
          label1 = labels$l1,
          label2 = labels$l2
        )
      )
    }
  )

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

  output$jaccard_heatmap <- renderPlot({
    dat <- loaded_jaccard_data()

    label1 <- if (!is.null(input$jaccard_label1) && nzchar(input$jaccard_label1)) input$jaccard_label1 else dat$l1
    label2 <- if (!is.null(input$jaccard_label2) && nzchar(input$jaccard_label2)) input$jaccard_label2 else dat$l2

    jaccard_heatmap_plot(
      meta1 = dat$meta1,
      meta2 = dat$meta2,
      name1 = label1,
      name2 = label2,
      threshold = input$jaccard_threshold
    )
  }, res = 110)

  output$download_jaccard_pdf <- downloadHandler(
    filename = function() {
      dat <- loaded_jaccard_data()
      paste0("cluster-similarity-heatmap-", dat$l1, "-vs-", dat$l2, ".pdf")
    },
    content = function(file) {
      dat <- loaded_jaccard_data()

      label1 <- if (!is.null(input$jaccard_label1) && nzchar(input$jaccard_label1)) input$jaccard_label1 else dat$l1
      label2 <- if (!is.null(input$jaccard_label2) && nzchar(input$jaccard_label2)) input$jaccard_label2 else dat$l2

      pdf_width <- 6.3
      pdf_height <- 6.3

      grDevices::pdf(file, width = pdf_width, height = pdf_height, onefile = TRUE)
      on.exit(grDevices::dev.off(), add = TRUE)
      print(jaccard_heatmap_plot(
        meta1 = dat$meta1,
        meta2 = dat$meta2,
        name1 = label1,
        name2 = label2,
        threshold = input$jaccard_threshold,
        for_pdf = TRUE
      ))
    }
  )

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

      pdf_width <- 7
      pdf_height <- 6

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

  output$jaccard_status <- renderText({
    dat <- loaded_jaccard_data()

    label1 <- if (!is.null(input$jaccard_label1) && nzchar(input$jaccard_label1)) input$jaccard_label1 else dat$l1
    label2 <- if (!is.null(input$jaccard_label2) && nzchar(input$jaccard_label2)) input$jaccard_label2 else dat$l2

    paste0(
      "Comparing ",
      dat$l1,
      " vs ",
      dat$l2,
      ". Axis labels: ",
      label1,
      " / ",
      label2,
      ". Heatmap cells at or above ",
      format(input$jaccard_threshold, trim = TRUE),
      " are annotated."
    )
  })
}

shinyApp(ui = ui, server = server)
