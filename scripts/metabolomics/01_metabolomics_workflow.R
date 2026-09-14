# =============================================================================
# Script:  01_metabolomics_workflow.R
# Purpose: Canonical preprocessing and QC/exploratory workflow for the
#          untargeted LC-MS/MS metabolomics feature tables (NEG, POS, and
#          COMB ESI modes) - missing-value and IQR filtering, median
#          normalization, log10 transformation, Pareto scaling, sample
#          distance heatmaps (Euclidean/Pearson/Spearman), PCA, t-SNE, UMAP,
#          per-time-point PCA/t-SNE/UMAP for COMB, and differential
#          accumulated metabolite (DAM) analysis with limma (see NOTE below).
#
#          Style notes kept from the original script:
#          - Colors follow the C. melo RNA-seq script's palette convention
#            (Mock/Inoculated from viridis; time points as sequential
#            viridis; heatmaps in viridis).
#          - DAM analysis: limma contrasts on median-normalized,
#            log10-transformed values (before Pareto scaling); threshold
#            adjusted p-value < 0.05 and |log2FC| > 1; DEG-style lollipop
#            summaries for (1) Inoculated vs Mock at each time point and
#            (2) sequential Inoculated comparisons.
#
#          NOTE: this limma-based DAM step differs from the one-way ANOVA
#          with Benjamini-Hochberg FDR correction described in the article's
#          Methods section for metabolomics. This has not been reconciled -
#          see the "Repository status" section of the README.
#
# Input:   data/metabolomics/002_METABO_NEG_AVG.txt
#          data/metabolomics/002_METABO_POS_AVG.txt
#          data/metabolomics/003_METABO_COMB_AVG.txt
# Output:  results/metabolomics/workflow/figures/
#          results/metabolomics/workflow/summaries/
#          results/metabolomics/workflow/metadata/
# Usage:   Rscript scripts/metabolomics/01_metabolomics_workflow.R
#          (run from the project root; see data/README.md to obtain inputs)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(readxl)
  library(Rtsne)
  library(uwot)
  library(viridis)
  library(grid)
  library(pheatmap)
  library(dendextend)
  library(limma)
})

############################################################
### Helper functions
############################################################

make_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

choose_tsne_perplexity <- function(sample_count, target_perplexity = 10) {
  max_valid <- floor((sample_count - 1) / 3)
  if (max_valid < 1) {
    return(NA_real_)
  }
  min(target_perplexity, max_valid)
}

get_canonical_tsne_settings <- function(sample_count) {
  list(
    perplexity = choose_tsne_perplexity(sample_count, target_perplexity = 10),
    theta = 0.5,
    max_iter = 2000
  )
}

compute_condition_score <- function(matrix, sampleinfo_metabo) {
  mock_idx <- sampleinfo_metabo$condition_long == "Mock"
  inoc_idx <- sampleinfo_metabo$condition_long == "Inoculated"
  
  if (sum(mock_idx) == 0 || sum(inoc_idx) == 0) {
    return(NA_real_)
  }
  
  mock_mat <- matrix[mock_idx, , drop = FALSE]
  inoc_mat <- matrix[inoc_idx, , drop = FALSE]
  
  mock_center <- colMeans(mock_mat)
  inoc_center <- colMeans(inoc_mat)
  between <- sqrt(sum((mock_center - inoc_center)^2))
  
  mock_within <- mean(sqrt(rowSums(sweep(mock_mat, 2, mock_center, "-")^2)))
  inoc_within <- mean(sqrt(rowSums(sweep(inoc_mat, 2, inoc_center, "-")^2)))
  within_mean <- mean(c(mock_within, inoc_within))
  
  if (!is.finite(within_mean) || within_mean == 0) {
    return(NA_real_)
  }
  
  between / within_mean
}

save_heatmap <- function(matrix, distance_obj, palette, filename, sample_annotation, annotation_colors, width = 10, height = 10) {
  pheatmap(
    matrix,
    clustering_distance_rows = distance_obj,
    clustering_distance_cols = distance_obj,
    annotation_col = sample_annotation,
    annotation_row = sample_annotation,
    annotation_colors = annotation_colors,
    color = palette,
    fontsize_col = 7,
    fontsize_row = 7,
    filename = filename,
    width = width,
    height = height
  )
}

############################################################
### Publication-style theme reused across metabolomics plots
############################################################

theme_metabo_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(size = 12, face = "bold", color = "black"),
      axis.text = element_text(size = 10, color = "black"),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.35, "cm"),
      plot.title = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
}

plot_pca <- function(
    pca_data,
    percent_var,
    color_mapping,
    shape_mapping,
    color_title,
    shape_title,
    filename,
    width = 8,
    height = 6,
    dpi = 300,
    additional_outputs = list()
) {
  p <- ggplot(
    pca_data,
    aes(
      x = PC1,
      y = PC2,
      color = .data[[color_mapping$column]],
      shape = .data[[shape_mapping$column]],
      label = sample
    )
  ) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(
      size = 3,
      max.overlaps = 20,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = "grey70",
      segment.size = 0.3,
      show.legend = FALSE
    ) +
    scale_color_manual(name = color_title, values = color_mapping$values) +
    scale_shape_manual(name = shape_title, values = shape_mapping$values) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(
      x = paste0("PC1 (", percent_var[1], "%)"),
      y = paste0("PC2 (", percent_var[2], "%)")
    ) +
    theme_metabo_publication()
  
  ggsave(
    filename = filename,
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    limitsize = FALSE
  )

  if (length(additional_outputs) > 0) {
    for (output in additional_outputs) {
      ggsave(
        filename = output$filename,
        plot = p,
        width = if (!is.null(output$width)) output$width else width,
        height = if (!is.null(output$height)) output$height else height,
        units = if (!is.null(output$units)) output$units else "in",
        dpi = if (!is.null(output$dpi)) output$dpi else dpi,
        limitsize = FALSE
      )
    }
  }
}

plot_tsne <- function(tsne_data, color_mapping, shape_mapping, color_title, shape_title, filename) {
  p <- ggplot(
    tsne_data,
    aes(
      x = tSNE1,
      y = tSNE2,
      color = .data[[color_mapping$column]],
      shape = .data[[shape_mapping$column]],
      label = sample
    )
  ) +
    geom_point(size = 3, alpha = 0.95) +
    geom_text_repel(
      size = 3,
      max.overlaps = 20,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = "grey70",
      segment.size = 0.3,
      show.legend = FALSE
    ) +
    scale_color_manual(name = color_title, values = color_mapping$values) +
    scale_shape_manual(name = shape_title, values = shape_mapping$values) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(
      x = "t-SNE 1",
      y = "t-SNE 2"
    ) +
    theme_metabo_publication()
  
  ggsave(filename, p, width = 8, height = 6, dpi = 300)
}

plot_umap <- function(umap_data, color_mapping, shape_mapping, color_title, shape_title, filename) {
  p <- ggplot(
    umap_data,
    aes(
      x = UMAP1,
      y = UMAP2,
      color = .data[[color_mapping$column]],
      shape = .data[[shape_mapping$column]],
      label = sample
    )
  ) +
    geom_point(size = 3, alpha = 0.95) +
    geom_text_repel(
      size = 3,
      max.overlaps = 20,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = "grey70",
      segment.size = 0.3,
      show.legend = FALSE
    ) +
    scale_color_manual(name = color_title, values = color_mapping$values) +
    scale_shape_manual(name = shape_title, values = shape_mapping$values) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_metabo_publication()
  
  ggsave(filename, p, width = 8, height = 6, dpi = 300)
}

run_umap <- function(matrix, seed = 123, target_neighbors = 10, min_dist = 0.3) {
  sample_count <- nrow(matrix)
  if (sample_count < 3) {
    return(NULL)
  }
  
  n_neighbors <- min(max(2, target_neighbors), sample_count - 1)
  
  set.seed(seed)
  umap_res <- uwot::umap(
    X = matrix,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    metric = "euclidean",
    verbose = FALSE,
    ret_model = FALSE
  )
  
  colnames(umap_res) <- c("UMAP1", "UMAP2")
  umap_res
}

############################################################
### Data reading
############################################################

read_metabo_dataset <- function(input_file, time_levels) {
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file)
  }

  file_ext <- tolower(tools::file_ext(input_file))

  if (file_ext %in% c("xlsx", "xlsm", "xls")) {
    metabo_raw <- readxl::read_excel(input_file, sheet = 1, .name_repair = "minimal")
    metabo_raw <- as.data.frame(metabo_raw, check.names = FALSE, stringsAsFactors = FALSE)

    if (nrow(metabo_raw) < 2 || ncol(metabo_raw) < 4) {
      stop("Input workbook does not contain the expected metabolomics matrix structure: ", input_file)
    }

    sample_header_raw <- colnames(metabo_raw)
    sample_header_clean <- sub("/.*$", "", sample_header_raw)
    sample_col_idx <- which(grepl("^\\d+h_[MI]_R[123]$", sample_header_clean))

    if (length(sample_col_idx) == 0) {
      stop("No sample columns were detected in workbook: ", input_file)
    }

    feature_names <- as.character(metabo_raw[[1]])
    feature_names[is.na(feature_names) | feature_names == ""] <- paste0("feature_", which(is.na(feature_names) | feature_names == ""))
    feature_names <- make.unique(feature_names)

    sample_names <- sample_header_clean[sample_col_idx]
    sample_labels <- sample_names

    numeric_feature_matrix <- metabo_raw[, sample_col_idx, drop = FALSE]
    numeric_feature_matrix <- as.data.frame(lapply(numeric_feature_matrix, as.numeric), check.names = FALSE)
    rownames(numeric_feature_matrix) <- feature_names

    sampleinfo_metabo <- data.frame(
      sample = sample_names,
      label = sample_labels,
      stringsAsFactors = FALSE
    ) %>%
      tidyr::separate(sample, into = c("time_point", "condition", "replicate"), sep = "_", remove = FALSE)

    if (anyNA(sampleinfo_metabo$time_point) || anyNA(sampleinfo_metabo$condition) || anyNA(sampleinfo_metabo$replicate)) {
      stop("Sample names must follow the pattern <time>_<condition>_<replicate> in workbook: ", input_file)
    }

    sampleinfo_metabo$time_point <- factor(sampleinfo_metabo$time_point, levels = time_levels)
    sampleinfo_metabo$condition <- factor(sampleinfo_metabo$condition, levels = c("M", "I"))
    sampleinfo_metabo$replicate <- factor(sampleinfo_metabo$replicate, levels = c("R1", "R2", "R3"))
    sampleinfo_metabo$condition_long <- factor(
      sampleinfo_metabo$condition,
      levels = c("M", "I"),
      labels = c("Mock", "Inoculated")
    )

    if (any(is.na(sampleinfo_metabo$time_point)) || any(is.na(sampleinfo_metabo$condition_long))) {
      stop("Unexpected time points or conditions found in workbook: ", input_file)
    }

    rownames(sampleinfo_metabo) <- sampleinfo_metabo$sample

    sample_annotation <- sampleinfo_metabo[, c("condition_long", "time_point"), drop = FALSE]
    colnames(sample_annotation) <- c("condition", "time_point")

    sample_feature_matrix_raw <- t(as.matrix(numeric_feature_matrix))
    rownames(sample_feature_matrix_raw) <- sample_names

    sample_feature_matrix_raw[!is.finite(sample_feature_matrix_raw)] <- NA_real_
    sample_feature_matrix_raw[sample_feature_matrix_raw < 0] <- NA_real_

    if (!identical(rownames(sample_feature_matrix_raw), rownames(sample_annotation))) {
      sample_annotation <- sample_annotation[rownames(sample_feature_matrix_raw), , drop = FALSE]
    }

    sampleinfo_metabo <- sampleinfo_metabo[rownames(sample_feature_matrix_raw), , drop = FALSE]
  } else {
    metabo_raw <- read.csv(
      input_file,
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    if (nrow(metabo_raw) < 2 || ncol(metabo_raw) < 2) {
      stop("Input file does not contain the expected metabolomics matrix structure: ", input_file)
    }

    if (!identical(as.character(metabo_raw[1, 1]), "Label")) {
      stop("The first data row must start with 'Label' in file: ", input_file)
    }

    feature_names <- make.unique(as.character(metabo_raw[-1, 1]))
    sample_names <- colnames(metabo_raw)[-1]
    sample_labels <- as.character(unlist(metabo_raw[1, -1]))

    numeric_feature_matrix <- metabo_raw[-1, -1, drop = FALSE]
    numeric_feature_matrix <- as.data.frame(lapply(numeric_feature_matrix, as.numeric), check.names = FALSE)
    rownames(numeric_feature_matrix) <- feature_names

    sampleinfo_metabo <- data.frame(
      sample = sample_names,
      label = sample_labels,
      stringsAsFactors = FALSE
    ) %>%
      tidyr::separate(sample, into = c("time_point", "condition", "replicate"), sep = "_", remove = FALSE)

    if (anyNA(sampleinfo_metabo$time_point) || anyNA(sampleinfo_metabo$condition) || anyNA(sampleinfo_metabo$replicate)) {
      stop("Sample names must follow the pattern <time>_<condition>_<replicate> in file: ", input_file)
    }

    sampleinfo_metabo$time_point <- factor(sampleinfo_metabo$time_point, levels = time_levels)
    sampleinfo_metabo$condition <- factor(sampleinfo_metabo$condition, levels = c("M", "I"))
    sampleinfo_metabo$replicate <- factor(sampleinfo_metabo$replicate, levels = c("R1", "R2", "R3"))
    sampleinfo_metabo$condition_long <- factor(
      sampleinfo_metabo$condition,
      levels = c("M", "I"),
      labels = c("Mock", "Inoculated")
    )

    if (any(is.na(sampleinfo_metabo$time_point)) || any(is.na(sampleinfo_metabo$condition_long))) {
      stop("Unexpected time points or conditions found in file: ", input_file)
    }

    rownames(sampleinfo_metabo) <- sampleinfo_metabo$sample

    sample_annotation <- sampleinfo_metabo[, c("condition_long", "time_point"), drop = FALSE]
    colnames(sample_annotation) <- c("condition", "time_point")

    sample_feature_matrix_raw <- t(as.matrix(numeric_feature_matrix))
    rownames(sample_feature_matrix_raw) <- sub("^X", "", rownames(sample_feature_matrix_raw))

    sample_feature_matrix_raw[!is.finite(sample_feature_matrix_raw)] <- NA_real_
    sample_feature_matrix_raw[sample_feature_matrix_raw < 0] <- NA_real_

    if (!identical(rownames(sample_feature_matrix_raw), rownames(sample_annotation))) {
      sample_annotation <- sample_annotation[rownames(sample_feature_matrix_raw), , drop = FALSE]
    }

    sampleinfo_metabo <- sampleinfo_metabo[rownames(sample_feature_matrix_raw), , drop = FALSE]
  }

  list(
    input_file = input_file,
    sample_feature_matrix_raw = sample_feature_matrix_raw,
    sampleinfo_metabo = sampleinfo_metabo,
    sample_annotation = sample_annotation,
    raw_feature_count = ncol(sample_feature_matrix_raw),
    sample_count = nrow(sample_feature_matrix_raw)
  )
}

############################################################
### Preprocessing functions
############################################################

filter_features_by_missing <- function(sample_feature_matrix_raw, missing_threshold = 0.50) {
  feature_missing_fraction <- colMeans(is.na(sample_feature_matrix_raw) | !is.finite(sample_feature_matrix_raw))
  keep_features <- feature_missing_fraction <= missing_threshold
  
  if (!any(keep_features)) {
    stop("No features remain after missing-value filtering.")
  }
  
  matrix_missing_filtered <- sample_feature_matrix_raw[, keep_features, drop = FALSE]
  
  list(
    matrix = matrix_missing_filtered,
    kept_count = ncol(matrix_missing_filtered),
    removed_count = sum(!keep_features)
  )
}

filter_features_by_iqr <- function(matrix_missing_filtered, fraction_to_remove = 0.40) {
  feature_iqr <- apply(matrix_missing_filtered, 2, IQR, na.rm = TRUE)
  feature_iqr[!is.finite(feature_iqr)] <- -Inf
  
  keep_fraction <- 1 - fraction_to_remove
  n_keep <- max(1, ceiling(length(feature_iqr) * keep_fraction))
  keep_idx <- order(feature_iqr, decreasing = TRUE)[seq_len(n_keep)]
  keep_idx <- sort(keep_idx)
  
  matrix_iqr_filtered <- matrix_missing_filtered[, keep_idx, drop = FALSE]
  
  list(
    matrix = matrix_iqr_filtered,
    kept_count = ncol(matrix_iqr_filtered),
    removed_count = length(feature_iqr) - ncol(matrix_iqr_filtered)
  )
}

median_normalize_samples <- function(matrix_iqr_filtered) {
  sample_medians <- apply(matrix_iqr_filtered, 1, function(x) median(x[is.finite(x)], na.rm = TRUE))
  invalid_medians <- !is.finite(sample_medians) | sample_medians <= 0
  
  if (any(invalid_medians)) {
    stop(
      "Median normalization failed because some samples have non-positive or undefined medians: ",
      paste(names(sample_medians)[invalid_medians], collapse = ", ")
    )
  }
  
  global_median <- median(sample_medians)
  matrix_median_normalized <- sweep(matrix_iqr_filtered, 1, sample_medians, "/")
  matrix_median_normalized <- matrix_median_normalized * global_median
  
  list(
    matrix = matrix_median_normalized,
    sample_medians = sample_medians,
    global_median = global_median
  )
}

log10_transform_matrix <- function(matrix_median_normalized, pseudocount = 1) {
  matrix_for_log10 <- matrix_median_normalized
  matrix_for_log10[matrix_for_log10 < 0] <- NA_real_
  
  matrix_log10_transformed <- log10(matrix_for_log10 + pseudocount)
  matrix_log10_transformed[!is.finite(matrix_log10_transformed)] <- NA_real_
  
  list(
    matrix = matrix_log10_transformed,
    pseudocount = pseudocount
  )
}

pareto_scale_features <- function(matrix_log10_transformed) {
  feature_means <- colMeans(matrix_log10_transformed, na.rm = TRUE)
  feature_sds <- apply(matrix_log10_transformed, 2, sd, na.rm = TRUE)
  valid_features <- is.finite(feature_means) & is.finite(feature_sds) & feature_sds > 0
  
  if (!any(valid_features)) {
    stop("No variable features remain for Pareto scaling.")
  }
  
  matrix_centered <- sweep(matrix_log10_transformed[, valid_features, drop = FALSE], 2, feature_means[valid_features], "-")
  pareto_denominators <- sqrt(feature_sds[valid_features])
  matrix_pareto_scaled <- sweep(matrix_centered, 2, pareto_denominators, "/")
  matrix_pareto_scaled[!is.finite(matrix_pareto_scaled)] <- NA_real_
  
  list(
    matrix = matrix_pareto_scaled,
    kept_count = ncol(matrix_pareto_scaled),
    removed_count = sum(!valid_features)
  )
}

finalize_processed_matrix <- function(matrix_pareto_scaled) {
  complete_features <- apply(matrix_pareto_scaled, 2, function(x) all(is.finite(x) & !is.na(x)))
  processed_matrix_final <- matrix_pareto_scaled[, complete_features, drop = FALSE]
  
  if (ncol(processed_matrix_final) == 0) {
    stop("No complete finite features remain after preprocessing.")
  }
  
  list(
    matrix = processed_matrix_final,
    kept_count = ncol(processed_matrix_final),
    removed_count = sum(!complete_features),
    complete_features = colnames(processed_matrix_final)
  )
}

preprocess_metabo_dataset <- function(raw_dataset, workflow_tag) {
  sample_feature_matrix_raw <- raw_dataset$sample_feature_matrix_raw
  
  missing_filter_result <- filter_features_by_missing(sample_feature_matrix_raw, missing_threshold = 0.50)
  matrix_missing_filtered <- missing_filter_result$matrix
  
  iqr_filter_result <- filter_features_by_iqr(matrix_missing_filtered, fraction_to_remove = 0.40)
  matrix_iqr_filtered <- iqr_filter_result$matrix
  
  median_normalization_result <- median_normalize_samples(matrix_iqr_filtered)
  matrix_median_normalized <- median_normalization_result$matrix
  
  log10_transform_result <- log10_transform_matrix(matrix_median_normalized, pseudocount = 1)
  matrix_log10_transformed <- log10_transform_result$matrix
  
  pareto_scaling_result <- pareto_scale_features(matrix_log10_transformed)
  matrix_pareto_scaled <- pareto_scaling_result$matrix
  
  final_cleanup_result <- finalize_processed_matrix(matrix_pareto_scaled)
  processed_matrix_final <- final_cleanup_result$matrix
  matrix_log10_final <- matrix_log10_transformed[
    rownames(processed_matrix_final),
    colnames(processed_matrix_final),
    drop = FALSE
  ]
  
  preprocessing_step_summary <- data.frame(
    workflow = workflow_tag,
    step_order = 1:7,
    step = c(
      "raw_input",
      "missing_filter_gt_50pct_removed",
      "iqr_filter_lowest_40pct_removed",
      "median_normalization",
      "log10_transformation",
      "pareto_scaling",
      "final_complete_finite_cleanup"
    ),
    feature_count = c(
      ncol(sample_feature_matrix_raw),
      ncol(matrix_missing_filtered),
      ncol(matrix_iqr_filtered),
      ncol(matrix_median_normalized),
      ncol(matrix_log10_transformed),
      ncol(matrix_pareto_scaled),
      ncol(processed_matrix_final)
    ),
    stringsAsFactors = FALSE
  )
  
  list(
    input_file = raw_dataset$input_file,
    processed_matrix_final = processed_matrix_final,
    matrix_log10_final = matrix_log10_final,
    sampleinfo_metabo = raw_dataset$sampleinfo_metabo,
    sample_annotation = raw_dataset$sample_annotation,
    raw_feature_count = ncol(sample_feature_matrix_raw),
    final_feature_count = ncol(processed_matrix_final),
    sample_count = nrow(processed_matrix_final),
    preprocessing_step_summary = preprocessing_step_summary
  )
}

############################################################
### DAM analysis: limma and DEG-style lollipop summaries
############################################################

run_limma_dam_contrast <- function(expr_matrix, group_vector, contrast_name) {
  group_vector <- factor(group_vector)
  design <- model.matrix(~ 0 + group_vector)
  colnames(design) <- levels(group_vector)
  
  fit <- limma::lmFit(t(expr_matrix), design)
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_name,
    levels = design
  )
  
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)
  
  res <- limma::topTable(
    fit2,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  
  res$metabolite <- rownames(res)
  res$log2FC <- res$logFC * log2(10)
  res$direction <- dplyr::case_when(
    !is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 & res$log2FC > 1 ~ "Up",
    !is.na(res$adj.P.Val) & res$adj.P.Val < 0.05 & res$log2FC < -1 ~ "Down",
    TRUE ~ "Not significant"
  )
  
  res[, c("metabolite", setdiff(colnames(res), "metabolite"))]
}

summarise_dam_result <- function(res, comparison_label, analysis_label) {
  sig_dams <- res %>%
    filter(
      !is.na(adj.P.Val),
      adj.P.Val < 0.05,
      abs(log2FC) > 1
    )
  
  data.frame(
    analysis = analysis_label,
    comparison = comparison_label,
    total = nrow(sig_dams),
    up = sum(sig_dams$log2FC > 1, na.rm = TRUE),
    down = sum(sig_dams$log2FC < -1, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

plot_dam_lollipop_summary <- function(dam_summary, analysis_label, comparison_order, x_label, filename) {
  plot_summary <- dam_summary %>%
    filter(analysis == analysis_label)
  
  if (nrow(plot_summary) == 0) {
    return(invisible(NULL))
  }
  
  dam_summary_long <- bind_rows(
    plot_summary %>%
      dplyr::select(comparison, up) %>%
      dplyr::mutate(direction = "Upregulated", count = up),
    plot_summary %>%
      dplyr::select(comparison, down) %>%
      dplyr::mutate(direction = "Downregulated", count = -down)
  ) %>%
    dplyr::select(comparison, direction, count)
  
  dam_summary_long$comparison_clean <- as.character(dam_summary_long$comparison)
  dam_summary_long$comparison_clean <- factor(
    dam_summary_long$comparison_clean,
    levels = comparison_order
  )
  
  dam_lollipop <- ggplot(
    dam_summary_long,
    aes(x = comparison_clean, y = count, color = direction)
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.5,
      color = "black"
    ) +
    geom_segment(
      aes(
        x = comparison_clean,
        xend = comparison_clean,
        y = 0,
        yend = count
      ),
      linewidth = 1,
      alpha = 0.8
    ) +
    geom_point(size = 4) +
    coord_flip() +
    scale_y_continuous(labels = abs) +
    scale_color_manual(
      values = c(
        Upregulated = "#D55E00",
        Downregulated = "#0072B2"
      )
    ) +
    labs(
      x = x_label,
      y = "DAMs"
    ) +
    theme_metabo_publication() +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 9)
    )
  
  ggsave(
    filename,
    dam_lollipop,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  invisible(dam_lollipop)
}

analyze_dams_metabo <- function(processed_dataset, dataset_tag, workflow_tag, figures_dir, summaries_dir, time_levels) {
  expr_matrix <- processed_dataset$matrix_log10_final
  sampleinfo_metabo <- processed_dataset$sampleinfo_metabo
  
  dam_figures_dir <- make_dir(file.path(figures_dir, dataset_tag, workflow_tag, "DAMs"))
  dam_tables_dir <- make_dir(file.path(summaries_dir, "DAMs", dataset_tag))
  
  all_results <- list()
  dam_summary_rows <- list()
  
  ############################################################
  ### A) Sequential DAMs in inoculated samples
  ############################################################
  
  inoc_info <- sampleinfo_metabo[sampleinfo_metabo$condition_long == "Inoculated", , drop = FALSE]
  inoc_matrix <- expr_matrix[rownames(inoc_info), , drop = FALSE]
  available_times <- time_levels[time_levels %in% as.character(inoc_info$time_point)]
  
  sequential_order <- character(0)
  
  if (length(available_times) >= 2) {
    sequential_pairs <- data.frame(
      current_time = available_times[-1],
      previous_time = available_times[-length(available_times)],
      stringsAsFactors = FALSE
    )
    
    sequential_order <- paste0(
      sequential_pairs$current_time,
      "_vs_",
      sequential_pairs$previous_time
    )
    
    for (i in seq_len(nrow(sequential_pairs))) {
      current_time <- sequential_pairs$current_time[i]
      previous_time <- sequential_pairs$previous_time[i]
      comparison_label <- paste0(current_time, "_vs_", previous_time)
      
      keep <- inoc_info$time_point %in% c(previous_time, current_time)
      sub_info <- inoc_info[keep, , drop = FALSE]
      sub_matrix <- inoc_matrix[rownames(sub_info), , drop = FALSE]
      
      if (nrow(sub_matrix) < 4 || length(unique(sub_info$time_point)) < 2) {
        next
      }
      
      group_vector <- paste0("T", sub_info$time_point)
      contrast_name <- paste0("T", current_time, "-T", previous_time)
      
      res <- run_limma_dam_contrast(
        expr_matrix = sub_matrix,
        group_vector = group_vector,
        contrast_name = contrast_name
      )
      
      res$analysis <- "sequential_inoculated"
      res$comparison <- comparison_label
      res$dataset <- dataset_tag
      
      write.csv(
        res,
        file.path(
          dam_tables_dir,
          paste0("DAM_", dataset_tag, "_sequential_inoculated_", comparison_label, "_all_metabolites.csv")
        ),
        row.names = FALSE
      )
      
      sig_res <- res %>%
        filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(log2FC) > 1)
      
      write.csv(
        sig_res,
        file.path(
          dam_tables_dir,
          paste0("DAM_", dataset_tag, "_sequential_inoculated_", comparison_label, "_significant.csv")
        ),
        row.names = FALSE
      )
      
      dam_summary_rows[[length(dam_summary_rows) + 1]] <- summarise_dam_result(
        res = res,
        comparison_label = comparison_label,
        analysis_label = "sequential_inoculated"
      )
      
      all_results[[length(all_results) + 1]] <- res
    }
  }
  
  ############################################################
  ### B) Inoculated vs Mock at each time point
  ############################################################
  
  inoculated_vs_mock_order <- character(0)
  
  for (time_point in time_levels) {
    keep <- sampleinfo_metabo$time_point == time_point
    sub_info <- sampleinfo_metabo[keep, , drop = FALSE]
    sub_matrix <- expr_matrix[rownames(sub_info), , drop = FALSE]
    
    if (nrow(sub_matrix) < 4 || length(unique(sub_info$condition_long)) < 2) {
      next
    }
    
    comparison_label <- paste0(time_point, "_Inoculated_vs_Mock")
    inoculated_vs_mock_order <- c(inoculated_vs_mock_order, comparison_label)
    
    res <- run_limma_dam_contrast(
      expr_matrix = sub_matrix,
      group_vector = sub_info$condition_long,
      contrast_name = "Inoculated-Mock"
    )
    
    res$analysis <- "inoculated_vs_mock"
    res$comparison <- comparison_label
    res$dataset <- dataset_tag
    
    write.csv(
      res,
      file.path(
        dam_tables_dir,
        paste0("DAM_", dataset_tag, "_", comparison_label, "_all_metabolites.csv")
      ),
      row.names = FALSE
    )
    
    sig_res <- res %>%
      filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(log2FC) > 1)
    
    write.csv(
      sig_res,
      file.path(
        dam_tables_dir,
        paste0("DAM_", dataset_tag, "_", comparison_label, "_significant.csv")
      ),
      row.names = FALSE
    )
    
    dam_summary_rows[[length(dam_summary_rows) + 1]] <- summarise_dam_result(
      res = res,
      comparison_label = comparison_label,
      analysis_label = "inoculated_vs_mock"
    )
    
    all_results[[length(all_results) + 1]] <- res
  }
  
  if (length(dam_summary_rows) == 0) {
    return(invisible(NULL))
  }
  
  dam_summary <- bind_rows(dam_summary_rows)
  dam_results <- bind_rows(all_results)
  
  write.csv(
    dam_summary,
    file.path(dam_tables_dir, paste0("DAM_", dataset_tag, "_summary.csv")),
    row.names = FALSE
  )
  
  write.csv(
    dam_results,
    file.path(dam_tables_dir, paste0("DAM_", dataset_tag, "_all_results.csv")),
    row.names = FALSE
  )
  
  plot_dam_lollipop_summary(
    dam_summary = dam_summary,
    analysis_label = "inoculated_vs_mock",
    comparison_order = inoculated_vs_mock_order,
    x_label = "Comparison (Inoculated vs Mock)",
    filename = file.path(
      dam_figures_dir,
      paste0("DAM_", dataset_tag, "_lollipop_inoculated_vs_mock.png")
    )
  )
  
  plot_dam_lollipop_summary(
    dam_summary = dam_summary,
    analysis_label = "sequential_inoculated",
    comparison_order = sequential_order,
    x_label = "Comparison (Sequential)",
    filename = file.path(
      dam_figures_dir,
      paste0("DAM_", dataset_tag, "_lollipop_sequential_inoculated.png")
    )
  )
  
  invisible(list(
    results = dam_results,
    summary = dam_summary
  ))
}

############################################################
### Main analysis functions
############################################################

analyze_processed_dataset <- function(processed_dataset, dataset_tag, workflow_tag, figures_dir, time_colors, condition_colors, annotation_colors, heatmap_colors) {
  message(
    "Processing ", dataset_tag,
    " | workflow: ", workflow_tag,
    " | source: ", processed_dataset$input_file
  )
  
  workflow_dir <- make_dir(file.path(figures_dir, dataset_tag, workflow_tag))
  tsne_dir <- make_dir(file.path(workflow_dir, "tSNE"))
  umap_dir <- make_dir(file.path(workflow_dir, "UMAP"))
  
  processed_matrix_final <- processed_dataset$processed_matrix_final
  sampleinfo_metabo <- processed_dataset$sampleinfo_metabo
  sample_annotation <- processed_dataset$sample_annotation
  
  ############################################################
  ### Euclidean sample-distance heatmap and dendrogram
  ############################################################
  
  sample_dists_euclidean <- dist(processed_matrix_final, method = "euclidean")
  sample_dist_matrix_euclidean <- as.matrix(sample_dists_euclidean)
  rownames(sample_dist_matrix_euclidean) <- rownames(processed_matrix_final)
  colnames(sample_dist_matrix_euclidean) <- rownames(processed_matrix_final)
  
  save_heatmap(
    sample_dist_matrix_euclidean,
    sample_dists_euclidean,
    heatmap_colors,
    file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_sample_distance_heatmap_euclidean.png")),
    sample_annotation,
    annotation_colors
  )
  
  sample_hclust <- hclust(sample_dists_euclidean, method = "ward.D2")
  dend <- as.dendrogram(sample_hclust)
  
  labels_colors(dend) <- condition_colors[
    as.character(sampleinfo_metabo$condition_long)
  ][order.dendrogram(dend)]
  
  dend <- set(dend, "labels_cex", 0.7)
  max_h <- attr(dend, "height")
  
  png(
    file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_sample_dendrogram_condition.png")),
    width = 2400,
    height = 3200,
    res = 300
  )
  
  par(mar = c(5, 12, 4, 4), xpd = TRUE)
  
  plot(
    dend,
    horiz = TRUE,
    main = "",
    xlab = "Height",
    xlim = c(max_h * 1.25, 0)
  )
  
  dev.off()
  
  ############################################################
  ### Pearson and Spearman sample-distance heatmaps
  ############################################################
  
  sample_cor_pearson <- cor(t(processed_matrix_final), method = "pearson")
  sample_dists_pearson <- as.dist(1 - sample_cor_pearson)
  sample_dist_matrix_pearson <- as.matrix(sample_dists_pearson)
  rownames(sample_dist_matrix_pearson) <- rownames(processed_matrix_final)
  colnames(sample_dist_matrix_pearson) <- rownames(processed_matrix_final)
  
  save_heatmap(
    sample_dist_matrix_pearson,
    sample_dists_pearson,
    heatmap_colors,
    file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_sample_distance_heatmap_pearson.png")),
    sample_annotation,
    annotation_colors
  )
  
  sample_cor_spearman <- cor(t(processed_matrix_final), method = "spearman")
  sample_dists_spearman <- as.dist(1 - sample_cor_spearman)
  sample_dist_matrix_spearman <- as.matrix(sample_dists_spearman)
  rownames(sample_dist_matrix_spearman) <- rownames(processed_matrix_final)
  colnames(sample_dist_matrix_spearman) <- rownames(processed_matrix_final)
  
  save_heatmap(
    sample_dist_matrix_spearman,
    sample_dists_spearman,
    heatmap_colors,
    file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_sample_distance_heatmap_spearman.png")),
    sample_annotation,
    annotation_colors
  )
  
  ############################################################
  ### PCA
  ############################################################
  
  pca_res <- prcomp(processed_matrix_final, center = FALSE, scale. = FALSE)
  pca_data <- as.data.frame(pca_res$x[, 1:2, drop = FALSE])
  pca_data$sample <- rownames(pca_data)
  pca_data <- left_join(pca_data, sampleinfo_metabo, by = "sample")
  
  percent_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))[1:2])
  
  plot_pca(
    pca_data = pca_data,
    percent_var = percent_var,
    color_mapping = list(column = "time_point", values = time_colors),
    shape_mapping = list(column = "condition_long", values = c(Mock = 17, Inoculated = 16)),
    color_title = "Time point",
    shape_title = "Condition",
    filename = file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_PCA.png")),
    width = 8,
    height = 6,
    additional_outputs = list(
      list(
        filename = file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_PCA.pdf")),
        width = 8,
        height = 6
      )
    )
  )
  
  plot_pca(
    pca_data = pca_data,
    percent_var = percent_var,
    color_mapping = list(column = "condition_long", values = condition_colors),
    shape_mapping = list(column = "condition_long", values = c(Mock = 17, Inoculated = 16)),
    color_title = "Condition",
    shape_title = "Condition",
    filename = file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_PCA_condition.png")),
    # PCA by condition: main horizontal version.
    width = 12,
    height = 5,
    additional_outputs = list(
      list(
        filename = file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_PCA_condition.pdf")),
        width = 12,
        height = 5
      ),
      list(
        filename = file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_PCA_condition.svg")),
        width = 12,
        height = 5
      ),
      list(
        filename = file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_PCA_condition_vertical.png")),
        width = 5,
        height = 12
      ),
      list(
        filename = file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_PCA_condition_square.png")),
        width = 8,
        height = 8
      )
    )
  )
  
  ############################################################
  ### t-SNE final canonical run
  ############################################################
  
  tsne_settings <- get_canonical_tsne_settings(nrow(processed_matrix_final))
  if (is.na(tsne_settings$perplexity)) {
    stop("The selected t-SNE perplexity is not valid for the number of samples in dataset: ", dataset_tag)
  }
  
  set.seed(123)
  tsne_res <- Rtsne(
    X = processed_matrix_final,
    dims = 2,
    perplexity = tsne_settings$perplexity,
    max_iter = tsne_settings$max_iter,
    theta = tsne_settings$theta,
    pca = TRUE,
    check_duplicates = FALSE
  )
  
  tsne_data <- as.data.frame(tsne_res$Y)
  colnames(tsne_data) <- c("tSNE1", "tSNE2")
  tsne_data$sample <- rownames(processed_matrix_final)
  tsne_data <- left_join(tsne_data, sampleinfo_metabo, by = "sample")
  
  write.csv(
    tsne_data,
    file.path(tsne_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_tSNE_final.csv")),
    row.names = FALSE
  )
  
  plot_tsne(
    tsne_data = tsne_data,
    color_mapping = list(column = "condition_long", values = condition_colors),
    shape_mapping = list(column = "condition_long", values = c(Mock = 17, Inoculated = 16)),
    color_title = "Condition",
    shape_title = "Condition",
    filename = file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_tSNE_condition_final.png"))
  )

  ############################################################
  ### UMAP final canonical run
  ############################################################
  
  umap_res <- run_umap(processed_matrix_final, seed = 123, target_neighbors = 10, min_dist = 0.3)
  if (!is.null(umap_res)) {
    umap_data <- as.data.frame(umap_res)
    umap_data$sample <- rownames(processed_matrix_final)
    umap_data <- left_join(umap_data, sampleinfo_metabo, by = "sample")
    
    write.csv(
      umap_data,
      file.path(umap_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_UMAP_final.csv")),
      row.names = FALSE
    )
    
    plot_umap(
      umap_data = umap_data,
      color_mapping = list(column = "condition_long", values = condition_colors),
      shape_mapping = list(column = "condition_long", values = c(Mock = 17, Inoculated = 16)),
      color_title = "Condition",
      shape_title = "Condition",
      filename = file.path(workflow_dir, paste0("Metabo_", dataset_tag, "_", workflow_tag, "_UMAP_condition_final.png"))
    )
  }
  
  message(
    "Finished ", dataset_tag,
    " | workflow: ", workflow_tag,
    " | samples: ", processed_dataset$sample_count,
    " | raw features: ", processed_dataset$raw_feature_count,
    " | final features: ", processed_dataset$final_feature_count
  )
  
  data.frame(
    dataset = dataset_tag,
    workflow = workflow_tag,
    sample_count = processed_dataset$sample_count,
    raw_feature_count = processed_dataset$raw_feature_count,
    final_feature_count = processed_dataset$final_feature_count,
    tsne_perplexity = tsne_settings$perplexity,
    tsne_theta = tsne_settings$theta,
    tsne_max_iter = tsne_settings$max_iter,
    umap_n_neighbors = if (!is.null(umap_res)) min(max(2, 10), nrow(processed_matrix_final) - 1) else NA_integer_,
    umap_min_dist = if (!is.null(umap_res)) 0.3 else NA_real_,
    source_file = processed_dataset$input_file,
    stringsAsFactors = FALSE
  )
}

analyze_combined_by_time <- function(processed_dataset, workflow_tag, figures_dir, summaries_dir, time_levels, condition_colors) {
  processed_matrix_final <- processed_dataset$processed_matrix_final
  sampleinfo_metabo <- processed_dataset$sampleinfo_metabo
  
  by_time_dir <- make_dir(file.path(figures_dir, "COMB", workflow_tag, "by_time"))
  time_points <- as.character(unique(sampleinfo_metabo$time_point))
  time_points <- time_points[!is.na(time_points)]
  time_points <- time_points[time_points %in% time_levels]
  time_points <- time_levels[time_levels %in% time_points]
  
  summary_rows <- list()
  
  for (time_point in time_points) {
    keep <- sampleinfo_metabo$time_point == time_point
    sub_info <- sampleinfo_metabo[keep, , drop = FALSE]
    sub_matrix <- processed_matrix_final[rownames(sub_info), , drop = FALSE]
    
    if (nrow(sub_matrix) < 4 || length(unique(sub_info$condition_long)) < 2) {
      next
    }
    
    message("Processing COMB by time | ", time_point, " | workflow: ", workflow_tag)
    
    time_dir <- make_dir(file.path(by_time_dir, time_point))
    tsne_dir <- make_dir(file.path(time_dir, "tSNE"))
    umap_dir <- make_dir(file.path(time_dir, "UMAP"))
    
    pca_res <- prcomp(sub_matrix, center = FALSE, scale. = FALSE)
    pca_data <- as.data.frame(pca_res$x[, 1:2, drop = FALSE])
    pca_data$sample <- rownames(pca_data)
    pca_data <- left_join(pca_data, sub_info, by = "sample")
    
    percent_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))[1:2])
    
    plot_pca(
      pca_data = pca_data,
      percent_var = percent_var,
      color_mapping = list(column = "condition_long", values = condition_colors),
      shape_mapping = list(column = "replicate", values = c(R1 = 16, R2 = 17, R3 = 15)),
      color_title = "Condition",
      shape_title = "Replicate",
      filename = file.path(time_dir, paste0("Metabo_COMB_", workflow_tag, "_", time_point, "_PCA.png")),
      width = 6,
      height = 5,
      additional_outputs = list(
        list(
          filename = file.path(time_dir, paste0("Metabo_COMB_", workflow_tag, "_", time_point, "_PCA.pdf")),
          width = 6,
          height = 5
        ),
        list(
          filename = file.path(time_dir, paste0("Metabo_COMB_", workflow_tag, "_", time_point, "_PCA.svg")),
          width = 6,
          height = 5
        )
      )
    )
    
    tsne_settings <- get_canonical_tsne_settings(nrow(sub_matrix))
    if (!is.na(tsne_settings$perplexity)) {
      set.seed(123)
      tsne_res <- Rtsne(
        X = sub_matrix,
        dims = 2,
        perplexity = tsne_settings$perplexity,
        max_iter = tsne_settings$max_iter,
        theta = tsne_settings$theta,
        pca = TRUE,
        check_duplicates = FALSE
      )
      
      tsne_data <- as.data.frame(tsne_res$Y)
      colnames(tsne_data) <- c("tSNE1", "tSNE2")
      tsne_data$sample <- rownames(sub_matrix)
      tsne_data <- left_join(tsne_data, sub_info, by = "sample")
      
      write.csv(
        tsne_data,
        file.path(tsne_dir, paste0("Metabo_COMB_", workflow_tag, "_", time_point, "_tSNE.csv")),
        row.names = FALSE
      )
      
      plot_tsne(
        tsne_data = tsne_data,
        color_mapping = list(column = "condition_long", values = condition_colors),
        shape_mapping = list(column = "replicate", values = c(R1 = 16, R2 = 17, R3 = 15)),
        color_title = "Condition",
        shape_title = "Replicate",
        filename = file.path(time_dir, paste0("Metabo_COMB_", workflow_tag, "_", time_point, "_tSNE.png"))
      )
    }

    umap_res <- run_umap(sub_matrix, seed = 123, target_neighbors = 10, min_dist = 0.3)
    if (!is.null(umap_res)) {
      umap_data <- as.data.frame(umap_res)
      umap_data$sample <- rownames(sub_matrix)
      umap_data <- left_join(umap_data, sub_info, by = "sample")
      
      write.csv(
        umap_data,
        file.path(umap_dir, paste0("Metabo_COMB_", workflow_tag, "_", time_point, "_UMAP.csv")),
        row.names = FALSE
      )
      
      plot_umap(
        umap_data = umap_data,
        color_mapping = list(column = "condition_long", values = condition_colors),
        shape_mapping = list(column = "replicate", values = c(R1 = 16, R2 = 17, R3 = 15)),
        color_title = "Condition",
        shape_title = "Replicate",
        filename = file.path(time_dir, paste0("Metabo_COMB_", workflow_tag, "_", time_point, "_UMAP.png"))
      )
    }
    
    umap_n_neighbors <- if (!is.null(umap_res)) min(max(2, 10), nrow(sub_matrix) - 1) else NA_integer_
    
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      time_point = time_point,
      sample_count = nrow(sub_matrix),
      feature_count = ncol(sub_matrix),
      tsne_perplexity = tsne_settings$perplexity,
      tsne_theta = tsne_settings$theta,
      tsne_max_iter = tsne_settings$max_iter,
      umap_n_neighbors = umap_n_neighbors,
      umap_min_dist = if (!is.null(umap_res)) 0.3 else NA_real_,
      condition_score = compute_condition_score(sub_matrix, sub_info),
      pca_pc1_percent = percent_var[1],
      pca_pc2_percent = percent_var[2],
      workflow = workflow_tag,
      source_file = processed_dataset$input_file,
      stringsAsFactors = FALSE
    )
  }
  
  if (length(summary_rows) == 0) {
    return(invisible(NULL))
  }
  
  summary_df <- bind_rows(summary_rows)
  write.csv(summary_df, file.path(summaries_dir, "COMB_by_time_summary.csv"), row.names = FALSE)
  write.csv(summary_df, file.path(by_time_dir, "COMB_by_time_summary.csv"), row.names = FALSE)
  
  invisible(summary_df)
}

############################################################
### Metadata
############################################################

write_metadata <- function(metadata_dir, dataset_specs, workflow_tag) {
  input_manifest <- bind_rows(lapply(dataset_specs, function(spec) {
    data.frame(
      dataset = spec$tag,
      input_file = spec$file,
      stringsAsFactors = FALSE
    )
  }))
  write.csv(input_manifest, file.path(metadata_dir, "input_manifest.csv"), row.names = FALSE)
  
  workflow_manifest <- data.frame(
    workflow = workflow_tag,
    step_order = 1:5,
    step = c(
      "Remove features with >50% missing values",
      "Remove the lowest 40% of features by IQR",
      "Median normalize samples",
      "Apply log10 transformation",
      "Apply Pareto scaling"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(workflow_manifest, file.path(metadata_dir, "preprocessing_manifest.csv"), row.names = FALSE)
  
  tsne_manifest <- data.frame(
    workflow = workflow_tag,
    target_perplexity = 10,
    theta = 0.5,
    max_iter = 2000,
    stringsAsFactors = FALSE
  )
  write.csv(tsne_manifest, file.path(metadata_dir, "tsne_manifest.csv"), row.names = FALSE)

  umap_manifest <- data.frame(
    workflow = workflow_tag,
    target_n_neighbors = 10,
    min_dist = 0.3,
    metric = "euclidean",
    stringsAsFactors = FALSE
  )
  write.csv(umap_manifest, file.path(metadata_dir, "umap_manifest.csv"), row.names = FALSE)
  
  writeLines(
    c(
      paste("Run timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      paste("R version:", R.version.string),
      paste("Workflow tag:", workflow_tag),
      "Canonical t-SNE settings: perplexity target=10, theta=0.5, max_iter=2000",
      "Canonical UMAP settings: n_neighbors target=10, min_dist=0.3, metric=euclidean",
      "DAM analysis: limma on median-normalized log10-transformed values before Pareto scaling",
      "DAM threshold: adjusted p-value < 0.05 and absolute log2FC > 1",
      "DAM lollipops: DEG-style diverging count plots for Upregulated and Downregulated DAMs"
    ),
    con = file.path(metadata_dir, "run_info.txt")
  )
  capture.output(sessionInfo(), file = file.path(metadata_dir, "session_info.txt"))
}

############################################################
### Run workflow
############################################################

run_metabo_workflow <- function(project_dir) {
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

  data_processed_dir <- file.path(project_dir, "data", "metabolomics")
  results_dir <- make_dir(file.path(project_dir, "results", "metabolomics", "workflow"))
  figures_dir <- make_dir(file.path(results_dir, "figures"))
  summaries_dir <- make_dir(file.path(results_dir, "summaries"))
  metadata_dir <- make_dir(file.path(results_dir, "metadata"))
  
  workflow_tag <- "missing50_iqr60_median_log10_pareto"
  
  dataset_specs <- list(
    list(file = file.path(data_processed_dir, "002_METABO_NEG_AVG.txt"), tag = "NEG"),
    list(file = file.path(data_processed_dir, "002_METABO_POS_AVG.txt"), tag = "POS"),
    list(file = file.path(data_processed_dir, "003_METABO_COMB_AVG.txt"), tag = "COMB")
  )
  
  time_levels <- c("12h", "24h", "36h", "48h", "60h", "72h", "84h", "96h", "108h", "120h", "132h", "144h")
  
  time_colors <- setNames(
    viridis::viridis(length(time_levels), option = "D"),
    time_levels
  )
  
  cmelo_green <- viridis::viridis(10, option = "D")[8]
  cmelo_inoculated <- viridis::viridis(10, option = "D")[4]
  
  condition_colors <- c(
    Mock = cmelo_green,
    Inoculated = cmelo_inoculated
  )
  
  annotation_colors <- list(
    condition = condition_colors,
    time_point = time_colors
  )
  
  heatmap_colors <- viridis::viridis(255, option = "D")
  
  analysis_summary_rows <- list()
  preprocessing_summary_rows <- list()
  combined_by_time_summary <- NULL
  dam_summary_rows <- list()
  
  for (spec in dataset_specs) {
    raw_dataset <- read_metabo_dataset(spec$file, time_levels)
    processed_dataset <- preprocess_metabo_dataset(raw_dataset, workflow_tag)
    
    preprocessing_step_summary <- processed_dataset$preprocessing_step_summary
    preprocessing_step_summary$dataset <- spec$tag
    preprocessing_summary_rows[[length(preprocessing_summary_rows) + 1]] <- preprocessing_step_summary
    
    analysis_summary_rows[[length(analysis_summary_rows) + 1]] <- analyze_processed_dataset(
      processed_dataset = processed_dataset,
      dataset_tag = spec$tag,
      workflow_tag = workflow_tag,
      figures_dir = figures_dir,
      time_colors = time_colors,
      condition_colors = condition_colors,
      annotation_colors = annotation_colors,
      heatmap_colors = heatmap_colors
    )
    
    dam_result <- analyze_dams_metabo(
      processed_dataset = processed_dataset,
      dataset_tag = spec$tag,
      workflow_tag = workflow_tag,
      figures_dir = figures_dir,
      summaries_dir = summaries_dir,
      time_levels = time_levels
    )
    
    if (!is.null(dam_result)) {
      dam_summary_rows[[length(dam_summary_rows) + 1]] <- dam_result$summary
    }
    
    if (spec$tag == "COMB") {
      combined_by_time_summary <- analyze_combined_by_time(
        processed_dataset = processed_dataset,
        workflow_tag = workflow_tag,
        figures_dir = figures_dir,
        summaries_dir = summaries_dir,
        time_levels = time_levels,
        condition_colors = condition_colors
      )
    }
  }
  
  preprocessing_summary_df <- bind_rows(preprocessing_summary_rows) %>%
    select(dataset, workflow, step_order, step, feature_count)
  write.csv(preprocessing_summary_df, file.path(summaries_dir, "preprocessing_step_summary.csv"), row.names = FALSE)
  
  analysis_summary_df <- bind_rows(analysis_summary_rows)
  write.csv(analysis_summary_df, file.path(summaries_dir, "analysis_input_summary.csv"), row.names = FALSE)
  
  if (length(dam_summary_rows) > 0) {
    dam_summary_all <- bind_rows(dam_summary_rows)
    write.csv(
      dam_summary_all,
      file.path(summaries_dir, "DAMs_all_datasets_summary.csv"),
      row.names = FALSE
    )
  } else {
    dam_summary_all <- NULL
  }
  
  write_metadata(metadata_dir, dataset_specs, workflow_tag)
  
  invisible(list(
    preprocessing_step_summary = preprocessing_summary_df,
    analysis_input_summary = analysis_summary_df,
    dam_summary = dam_summary_all,
    comb_by_time_summary = combined_by_time_summary
  ))
}

############################################################
### Script entry point
### Assumes the working directory is the project root (see Usage in the
### header above), consistent with the other scripts in this repository.
############################################################

if (sys.nframe() == 0) {
  run_metabo_workflow(project_dir = ".")
}
