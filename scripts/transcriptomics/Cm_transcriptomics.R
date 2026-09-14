############################################################
### Cucumis melo RNA-seq time-course analysis
### Publication-oriented integrated script
###
### Includes:
### - TPM export
### - Sample distances, dendrogram, PCA, t-SNE, UMAP
### - DESeq2 with condition * time_point design
### - Inoculated vs Mock contrasts at each time point
### - Distribution of expression changes (violin plot)
### - Temporal clustering of DEGs (k-means, Inoculated trajectory)
### - Functional enrichment analysis (GO / KEGG)
###
### Color strategy:
### - Mock: C. melo green from viridis
### - Inoculated: darker blue-green from viridis
### - Up/down DEG summaries: accessible signed-effect colors from plasma
############################################################

############################################################
### 1. Load required libraries
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(viridis)
  library(ggrepel)
  library(Rtsne)
  library(uwot)
  library(dendextend)
  library(tidyr)
  library(grid)
})

has_clusterprofiler <- requireNamespace("clusterProfiler", quietly = TRUE)
has_enrichplot <- requireNamespace("enrichplot", quietly = TRUE)

############################################################
### 2. Define input files and output directories
############################################################

project_root <- if (dir.exists(file.path(getwd(), "data", "processed"))) {
  file.path(getwd(), "data", "processed")
} else {
  getwd()
}

counts_file <- file.path(project_root, "01_Transcriptomics_Cm_counts.txt")
sampleinfo_file <- file.path(project_root, "02_Transcriptomics_Cm_sample_info.txt")
annotation_file <- file.path(project_root, "Cm_functional_annotation.tsv")   # Optional

results_root <- file.path(project_root, "results_melon")
figures_dir <- file.path(results_root, "figures")
tables_dir <- file.path(results_root, "tables")
deseq_dir <- file.path(tables_dir, "deseq2")
tsne_dir <- file.path(tables_dir, "tsne")
umap_dir <- file.path(tables_dir, "umap")
contrasts_dir <- file.path(deseq_dir, "condition_by_time")
clustering_dir <- file.path(tables_dir, "clustering")
enrichment_dir <- file.path(tables_dir, "enrichment")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deseq_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tsne_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(contrasts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(clustering_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrichment_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
### 3. Read count matrix and sample metadata
############################################################

counts <- read.table(
  counts_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sampleinfo <- read.table(
  sampleinfo_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

############################################################
### 4. Detect gene identifier column and validate required columns
############################################################

if ("Geneid" %in% colnames(counts)) {
  gene_id_col <- "Geneid"
} else if ("Name" %in% colnames(counts)) {
  gene_id_col <- "Name"
} else {
  stop("No gene identifier column found. Expected either 'Geneid' or 'Name'.")
}

required_sampleinfo_cols <- c("sample", "time_point", "replicate", "condition")
missing_sampleinfo_cols <- setdiff(required_sampleinfo_cols, colnames(sampleinfo))

if (length(missing_sampleinfo_cols) > 0) {
  stop(
    paste(
      "The sample metadata file is missing required columns:",
      paste(missing_sampleinfo_cols, collapse = ", ")
    )
  )
}

if (!"Length" %in% colnames(counts)) {
  stop("The count matrix file is missing required column: Length")
}

############################################################
### 5. Validate sample names and reorder metadata
############################################################

count_sample_names <- colnames(counts)[7:ncol(counts)]

missing_in_sampleinfo <- setdiff(count_sample_names, sampleinfo$sample)
missing_in_counts <- setdiff(sampleinfo$sample, count_sample_names)

if (length(missing_in_sampleinfo) > 0 || length(missing_in_counts) > 0) {
  stop("Sample names do not match between the count matrix and sample metadata.")
}

sampleinfo <- sampleinfo[match(count_sample_names, sampleinfo$sample), ]
rownames(sampleinfo) <- sampleinfo$sample

############################################################
### 6. Format experimental variables
############################################################

time_numeric <- as.integer(sub("h$", "", sampleinfo$time_point))
time_levels <- paste0(sort(unique(time_numeric)), "h")

sampleinfo$time_point <- factor(sampleinfo$time_point, levels = time_levels)
sampleinfo$replicate <- factor(sampleinfo$replicate)
sampleinfo$condition <- factor(sampleinfo$condition, levels = c("Mock", "Inoculated"))

if ("infection_phase" %in% colnames(sampleinfo)) {
  sampleinfo$infection_phase <- factor(
    sampleinfo$infection_phase,
    levels = c("Early", "Intermediate", "Late")
  )
}

sampleinfo$group <- factor(
  paste(sampleinfo$condition, sampleinfo$time_point, sep = "_"),
  levels = c(
    paste0("Mock_", time_levels),
    paste0("Inoculated_", time_levels)
  )
)

############################################################
### 7. Build count matrix and validate ordering
############################################################

count_data <- as.matrix(counts[, count_sample_names])
rownames(count_data) <- counts[[gene_id_col]]
storage.mode(count_data) <- "integer"

if (!identical(colnames(count_data), rownames(sampleinfo))) {
  stop("Sample order is not identical between the count matrix and sample metadata.")
}

if (any(duplicated(rownames(count_data)))) {
  stop("Duplicated gene identifiers detected in the count matrix.")
}

############################################################
### 8. Extract gene lengths and validate TPM inputs
############################################################

gene_lengths <- counts$Length
names(gene_lengths) <- counts[[gene_id_col]]

if (any(is.na(gene_lengths)) || any(gene_lengths <= 0)) {
  stop("Invalid gene lengths were detected for TPM calculation.")
}

############################################################
### 9. Define annotations, colors, and reusable plot theme
############################################################

### Condition colors based on viridis.
### Mock keeps the C. melo green identity.
### Inoculated is distinct but still within viridis, avoiding red/green contrast.
cmelo_green <- viridis::viridis(10, option = "D")[8]
cmelo_inoculated <- viridis::viridis(10, option = "D")[4]

condition_colors <- c(
  Mock = cmelo_green,
  Inoculated = cmelo_inoculated
)

### Time points as sequential viridis colors.
time_colors <- setNames(
  viridis::viridis(length(time_levels), option = "D"),
  time_levels
)

annotation_colors <- list(
  condition = condition_colors,
  time_point = time_colors
)

sample_annotation <- sampleinfo[, c("time_point", "condition"), drop = FALSE]

### Sequential palette for distances / expression intensity.
heatmap_colors <- viridis::viridis(255, option = "D")

### Signed effect palette for DEG direction.
### Kept separate from condition colors because up/down is a different semantic variable.
deg_direction_colors <- c(
  Upregulated = "#D62728",
  Downregulated = "#1F77B4"
)
### Publication-style theme reused across plots.
theme_cmelo_publication <- function(base_size = 12) {
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

############################################################
### 10. Calculate TPM values and export table
############################################################

gene_lengths_kb <- gene_lengths / 1000
rpk <- sweep(count_data, 1, gene_lengths_kb, "/")
scaling_factors <- colSums(rpk)

if (any(scaling_factors == 0)) {
  stop("At least one sample has zero total RPK, TPM calculation cannot proceed.")
}

tpm <- sweep(rpk, 2, scaling_factors, "/") * 1e6

counts_tpm <- cbind(counts[, 1:6], as.data.frame(tpm, check.names = FALSE))

write.table(
  counts_tpm,
  file.path(tables_dir, "Cmelo_counts_TPM.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

############################################################
### 11. Build DESeq2 dataset
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sampleinfo,
  design = ~ condition * time_point
)

############################################################
### 12. Filter low-expression genes and normalize counts
############################################################

keep_genes <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep_genes, ]

dds <- estimateSizeFactors(dds)

message(sum(keep_genes), " genes retained after low-count filtering.")

############################################################
### 13. Apply variance-stabilizing transformation
###     and save/reuse transformed object
############################################################

rld_file <- file.path(tables_dir, "Cmelo_vst_rld_object.rds")

if (file.exists(rld_file)) {
  message("Loading existing vst-transformed object from: ", rld_file)
  rld <- readRDS(rld_file)
} else {
  message("Computing vst transformation...")
  rld <- vst(dds, blind = FALSE)
  saveRDS(rld, file = rld_file)
  message("vst-transformed object saved to: ", rld_file)
}

############################################################
### 14. Define reusable helper functions
############################################################

save_heatmap <- function(matrix, distance_obj, palette, filename, width = 10, height = 10) {
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

filter_and_save_deg <- function(
    res,
    filename,
    comparison_label,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
) {
  res_df <- as.data.frame(res)
  res_df$Geneid <- rownames(res_df)
  res_df <- res_df[, c("Geneid", setdiff(colnames(res_df), "Geneid"))]
  
  sig_deg <- res_df[
    !is.na(res_df$padj) &
      res_df$padj < padj_cutoff &
      abs(res_df$log2FoldChange) > lfc_cutoff,
  ]
  
  sig_up <- sig_deg[
    sig_deg$log2FoldChange > lfc_cutoff,
  ]
  
  sig_down <- sig_deg[
    sig_deg$log2FoldChange < -lfc_cutoff,
  ]
  
  write.csv(
    sig_deg,
    filename,
    row.names = FALSE
  )
  
  base_name <- tools::file_path_sans_ext(filename)
  
  write.table(
    sig_up$Geneid,
    file = paste0(base_name, "_UP.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  write.table(
    sig_down$Geneid,
    file = paste0(base_name, "_DOWN.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  write.csv(
    sig_up,
    file = paste0(base_name, "_UP.csv"),
    row.names = FALSE
  )
  
  write.csv(
    sig_down,
    file = paste0(base_name, "_DOWN.csv"),
    row.names = FALSE
  )
  
  data.frame(
    comparison = comparison_label,
    total = nrow(sig_deg),
    up = nrow(sig_up),
    down = nrow(sig_down),
    stringsAsFactors = FALSE
  )
}

run_term_enrichment <- function(
    gene_ids,
    universe_ids,
    term2gene,
    out_csv,
    out_png,
    plot_title,
    top_n = 15
) {
  if (!has_clusterprofiler || !has_enrichplot) {
    message("clusterProfiler and/or enrichplot not installed. Skipping enrichment.")
    return(NULL)
  }
  
  if (length(gene_ids) < 2) {
    return(NULL)
  }
  
  gene_ids <- unique(gene_ids)
  universe_ids <- unique(universe_ids)
  term2gene <- unique(term2gene)
  
  term2gene <- term2gene[complete.cases(term2gene), , drop = FALSE]
  colnames(term2gene) <- c("term", "gene")
  
  if (nrow(term2gene) == 0) {
    return(NULL)
  }
  
  enr <- clusterProfiler::enricher(
    gene = gene_ids,
    universe = universe_ids,
    TERM2GENE = term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
  
  if (is.null(enr)) {
    return(NULL)
  }
  
  enr_df <- as.data.frame(enr)
  if (nrow(enr_df) == 0) {
    return(NULL)
  }
  
  write.csv(enr_df, out_csv, row.names = FALSE)
  
  png(out_png, width = 2400, height = 1800, res = 300)
  print(
    enrichplot::dotplot(enr, showCategory = min(top_n, nrow(enr_df))) +
      ggtitle(plot_title)
  )
  dev.off()
  
  invisible(enr_df)
}

############################################################
### 15. Sample-to-sample distance heatmap and dendrogram
############################################################

sample_dists_euclidean <- dist(t(assay(rld)), method = "euclidean")
sample_dist_matrix_euclidean <- as.matrix(sample_dists_euclidean)
rownames(sample_dist_matrix_euclidean) <- colnames(rld)
colnames(sample_dist_matrix_euclidean) <- colnames(rld)

save_heatmap(
  sample_dist_matrix_euclidean,
  sample_dists_euclidean,
  heatmap_colors,
  file.path(figures_dir, "Cmelo_sample_distance_heatmap_euclidean.png")
)

sample_dists <- dist(t(assay(rld)), method = "euclidean")
hc <- hclust(sample_dists, method = "ward.D2")
dend <- as.dendrogram(hc)

labels_colors(dend) <- condition_colors[
  as.character(sampleinfo$condition)
][order.dendrogram(dend)]

dend <- set(dend, "labels_cex", 0.7)
max_h <- attr(dend, "height")

png(
  file.path(figures_dir, "Cmelo_sample_dendrogram_condition.png"),
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
### 16. PCA colored by condition
############################################################

pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)

pca_data$condition <- factor(pca_data$condition, levels = c("Mock", "Inoculated"))
pca_data$sample <- rownames(pca_data)

percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = condition, shape = condition, label = sample)
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
  scale_color_manual(
    name = "Condition",
    values = condition_colors
  ) +
  scale_shape_manual(
    name = "Condition",
    values = c(
      Mock = 17,
      Inoculated = 16
    )
  ) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_y_continuous(expand = expansion(mult = 0.02)) +
  labs(
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  ) +
  theme_cmelo_publication()

ggsave(
  file.path(figures_dir, "Cmelo_PCA_condition.png"),
  pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

############################################################
### 17. Run t-SNE for exploratory visualization
############################################################

tsne_input <- t(assay(rld))

if ((nrow(tsne_input) - 1) < (3 * 5)) {
  message("The selected t-SNE perplexity is not valid for the number of samples. Skipping t-SNE.")
} else {
  
  set.seed(123)
  
  tsne_res <- Rtsne(
    X = tsne_input,
    dims = 2,
    perplexity = 5,
    max_iter = 2000,
    theta = 0.5,
    pca = TRUE,
    check_duplicates = FALSE
  )
  
  tsne_data <- as.data.frame(tsne_res$Y)
  colnames(tsne_data) <- c("tSNE1", "tSNE2")
  tsne_data$sample <- colnames(rld)
  
  tsne_data$condition <- factor(
    sampleinfo[tsne_data$sample, "condition"],
    levels = c("Mock", "Inoculated")
  )
  
  tsne_data$time_point <- factor(
    sampleinfo[tsne_data$sample, "time_point"],
    levels = time_levels
  )
  
  tsne_data$replicate <- sampleinfo[tsne_data$sample, "replicate"]
  
  write.csv(
    tsne_data,
    file.path(tsne_dir, "Cmelo_tSNE_paper.csv"),
    row.names = FALSE
  )
  
  tsne_plot <- ggplot(
    tsne_data,
    aes(x = tSNE1, y = tSNE2, color = condition, shape = condition, label = sample)
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
    scale_color_manual(
      name = "Condition",
      values = condition_colors
    ) +
    scale_shape_manual(
      name = "Condition",
      values = c(
        Mock = 17,
        Inoculated = 16
      )
    ) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(
      x = "t-SNE 1",
      y = "t-SNE 2"
    ) +
    theme_cmelo_publication()
  
  ggsave(
    file.path(figures_dir, "Cmelo_tSNE_condition_labeled.png"),
    tsne_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

############################################################
### 18. Run UMAP for exploratory visualization
############################################################

umap_input <- t(assay(rld))
n_samples <- nrow(umap_input)
umap_n_neighbors <- min(15, n_samples - 1)
umap_n_neighbors <- max(2, umap_n_neighbors)

if (n_samples < 3) {
  message("Not enough samples for UMAP. Skipping UMAP.")
} else {
  
  set.seed(123)
  
  umap_res <- uwot::umap(
    X = umap_input,
    n_neighbors = umap_n_neighbors,
    min_dist = 0.3,
    metric = "euclidean",
    verbose = FALSE,
    ret_model = FALSE
  )
  
  umap_data <- as.data.frame(umap_res)
  colnames(umap_data) <- c("UMAP1", "UMAP2")
  umap_data$sample <- colnames(rld)
  
  umap_data$condition <- factor(
    sampleinfo[umap_data$sample, "condition"],
    levels = c("Mock", "Inoculated")
  )
  
  umap_data$time_point <- factor(
    sampleinfo[umap_data$sample, "time_point"],
    levels = time_levels
  )
  
  umap_data$replicate <- sampleinfo[umap_data$sample, "replicate"]
  
  write.csv(
    umap_data,
    file.path(umap_dir, "Cmelo_UMAP_paper.csv"),
    row.names = FALSE
  )
  
  umap_plot <- ggplot(
    umap_data,
    aes(x = UMAP1, y = UMAP2, color = condition, shape = condition, label = sample)
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
    scale_color_manual(
      name = "Condition",
      values = condition_colors
    ) +
    scale_shape_manual(
      name = "Condition",
      values = c(
        Mock = 17,
        Inoculated = 16
      )
    ) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_cmelo_publication()
  
  ggsave(
    file.path(figures_dir, "Cmelo_UMAP_condition_labeled.png"),
    umap_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

############################################################
### 19. Run DESeq2 differential expression analysis
############################################################

dds_file <- file.path(tables_dir, "Cmelo_dds_DESeq_object.rds")

if (file.exists(dds_file)) {
  message("Loading existing DESeq2 object from: ", dds_file)
  dds <- readRDS(dds_file)
} else {
  message("Running DESeq2...")
  dds <- DESeq(dds)
  
  saveRDS(
    dds,
    file = dds_file
  )
  
  message("DESeq2 object saved to: ", dds_file)
}

############################################################
### 19. Save model coefficient names
############################################################

write.table(
  resultsNames(dds),
  file = file.path(deseq_dir, "Cmelo_resultsNames.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

############################################################
### 20. Compare Inoculated versus Mock at each time point
### Robust implementation using a group-based design
############################################################

dds_group_file <- file.path(tables_dir, "Cmelo_dds_group_object.rds")

if (file.exists(dds_group_file)) {
  message("Loading existing group-based DESeq2 object from: ", dds_group_file)
  dds_group <- readRDS(dds_group_file)
} else {
  
  dds_group <- DESeqDataSetFromMatrix(
    countData = counts(dds, normalized = FALSE),
    colData = as.data.frame(colData(dds)),
    design = ~ group
  )
  
  dds_group <- DESeq(dds_group)
  
  saveRDS(
    dds_group,
    file = dds_group_file
  )
  
  message("Group-based DESeq2 object saved to: ", dds_group_file)
}

deg_summary_condition <- do.call(
  rbind,
  lapply(time_levels, function(tp) {
    
    res <- results(
      dds_group,
      contrast = c("group", paste0("Inoculated_", tp), paste0("Mock_", tp))
    )
    
    comparison_label <- paste0("Inoculated_vs_Mock_", tp)
    output_file <- file.path(contrasts_dir, paste0("DEG_", comparison_label, ".csv"))
    
    filter_and_save_deg(res, output_file, comparison_label)
  })
)

write.csv(
  deg_summary_condition,
  file.path(contrasts_dir, "DEG_summary_Inoculated_vs_Mock_by_time.csv"),
  row.names = FALSE
)

############################################################
### 21. Plot DEG summary for Inoculated vs Mock at each time
############################################################

deg_summary_condition_long <- bind_rows(
  deg_summary_condition %>%
    dplyr::select(comparison, up) %>%
    dplyr::mutate(direction = "Upregulated", count = up),
  deg_summary_condition %>%
    dplyr::select(comparison, down) %>%
    dplyr::mutate(direction = "Downregulated", count = -down)
) %>%
  dplyr::select(comparison, direction, count)

deg_summary_condition_long$comparison <- factor(
  deg_summary_condition_long$comparison,
  levels = paste0("Inoculated_vs_Mock_", time_levels)
)

deg_barplot <- ggplot(
  deg_summary_condition_long,
  aes(x = comparison, y = count, fill = direction)
) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = deg_direction_colors
  ) +
  labs(
    title = "Differentially expressed genes: Inoculated vs Mock at each time point",
    x = "Comparison",
    y = "Number of DEGs"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_blank()
  )

ggsave(
  file.path(figures_dir, "Cmelo_DEG_summary_Inoculated_vs_Mock.png"),
  deg_barplot,
  width = 10,
  height = 7,
  dpi = 300
)

############################################################
### 22. Distribution of expression changes
###     Inoculated vs Mock at each time point
############################################################

lfc_list_condition <- lapply(time_levels, function(tp) {
  
  res <- results(
    dds_group,
    contrast = c("group", paste0("Inoculated_", tp), paste0("Mock_", tp))
  )
  
  res_df <- as.data.frame(res)
  res_df$Geneid <- rownames(res_df)
  res_df$comparison <- paste0("Inoculated_vs_Mock_", tp)
  res_df
})

lfc_data_condition <- bind_rows(lfc_list_condition) %>%
  dplyr::filter(!is.na(log2FoldChange))

if (nrow(lfc_data_condition) == 0) {
  
  message("No valid log2FoldChange values available for melon violin plot.")
  
} else {
  
  comparison_order_condition <- paste0("Inoculated_vs_Mock_", time_levels)
  
  lfc_data_condition$comparison <- factor(
    lfc_data_condition$comparison,
    levels = comparison_order_condition
  )
  
  write.csv(
    lfc_data_condition,
    file.path(contrasts_dir, "LFC_Inoculated_vs_Mock_all_timepoints.csv"),
    row.names = FALSE
  )
  
  deg_counts_labels_condition <- lfc_data_condition %>%
    dplyr::mutate(
      is_deg = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1
    ) %>%
    dplyr::group_by(comparison) %>%
    dplyr::summarise(
      n_deg = sum(is_deg, na.rm = TRUE),
      .groups = "drop"
    )
  
  deg_counts_labels_condition$label <- paste0("n=", deg_counts_labels_condition$n_deg)
  
  lfc_violin_condition <- ggplot(
    lfc_data_condition,
    aes(x = comparison, y = log2FoldChange)
  ) +
    geom_violin(
      fill = cmelo_inoculated,
      color = "black",
      trim = TRUE,
      linewidth = 0.4,
      alpha = 0.65
    ) +
    geom_boxplot(
      width = 0.12,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.4
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey30",
      linewidth = 0.6
    ) +
    geom_text(
      data = deg_counts_labels_condition,
      aes(x = comparison, y = 3.6, label = label),
      inherit.aes = FALSE,
      size = 3.2,
      fontface = "bold"
    ) +
    coord_cartesian(ylim = c(-4, 4)) +
    labs(
      title = "Distribution of expression changes: Inoculated vs Mock",
      x = "Comparison",
      y = "log2 fold change"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  print(lfc_violin_condition)
  
  ggsave(
    file.path(figures_dir, "Cmelo_log2FC_violin_Inoculated_vs_Mock.png"),
    plot = lfc_violin_condition,
    width = 10,
    height = 6.5,
    dpi = 300
  )
}

############################################################
### 23. Temporal clustering of DEGs
###     K-means on inoculation-responsive genes
############################################################

deg_files_condition <- list.files(
  contrasts_dir,
  pattern = "^DEG_Inoculated_vs_Mock_[0-9]+h\\.csv$",
  full.names = TRUE
)

if (length(deg_files_condition) == 0) {
  
  message("No DEG files found for Inoculated vs Mock contrasts. Skipping clustering.")
  
} else {
  
  deg_lists_condition <- lapply(deg_files_condition, read.csv)
  
  deg_genes_condition <- unique(
    unlist(lapply(deg_lists_condition, function(x) x$Geneid))
  )
  
  deg_genes_condition <- intersect(deg_genes_condition, rownames(rld))
  
  message(length(deg_genes_condition), " unique DEGs available for temporal clustering.")
  
  if (length(deg_genes_condition) < 2) {
    
    message("Not enough dynamic genes for melon k-means clustering. Skipping clustering.")
    
  } else {
    
    rld_mat <- assay(rld)
    
    inoculated_samples <- rownames(sampleinfo)[sampleinfo$condition == "Inoculated"]
    
    if (length(inoculated_samples) < 2) {
      
      message("Not enough inoculated samples for temporal clustering. Skipping clustering.")
      
    } else {
      
      rld_deg_inoc <- rld_mat[deg_genes_condition, inoculated_samples, drop = FALSE]
      sampleinfo_inoc <- sampleinfo[inoculated_samples, , drop = FALSE]
      
      time_means_inoc <- sapply(levels(sampleinfo_inoc$time_point), function(tp) {
        rowMeans(rld_deg_inoc[, sampleinfo_inoc$time_point == tp, drop = FALSE])
      })
      
      time_means_inoc <- as.matrix(time_means_inoc)
      colnames(time_means_inoc) <- levels(sampleinfo_inoc$time_point)
      
      time_means_scaled <- t(scale(t(time_means_inoc)))
      time_means_scaled <- time_means_scaled[complete.cases(time_means_scaled), , drop = FALSE]
      
      n_genes_cluster <- nrow(time_means_scaled)
      k <- min(3, n_genes_cluster)
      
      if (k < 2) {
        
        message("Not enough genes remaining after scaling for melon k-means clustering. Skipping.")
        
      } else {
        
        set.seed(123)
        
        km <- kmeans(
          time_means_scaled,
          centers = k,
          nstart = 50
        )
        
        cluster_df <- as.data.frame(time_means_scaled)
        cluster_df$Geneid <- rownames(time_means_scaled)
        cluster_df$cluster <- factor(km$cluster)
        
        write.csv(
          cluster_df[, c("Geneid", "cluster")],
          file.path(clustering_dir, "Cmelo_temporal_gene_clusters.csv"),
          row.names = FALSE
        )
        
        cluster_centers <- as.data.frame(km$centers)
        cluster_centers$cluster <- factor(seq_len(nrow(cluster_centers)))
        
        write.csv(
          cluster_centers,
          file.path(clustering_dir, "Cmelo_temporal_cluster_centers.csv"),
          row.names = FALSE
        )
        
        cluster_long <- cluster_df %>%
          tidyr::pivot_longer(
            cols = all_of(levels(sampleinfo_inoc$time_point)),
            names_to = "time_point",
            values_to = "scaled_expression"
          )
        
        cluster_long$time_point <- factor(
          cluster_long$time_point,
          levels = levels(sampleinfo_inoc$time_point)
        )
        
        cluster_plot <- ggplot(
          cluster_long,
          aes(x = time_point, y = scaled_expression, group = Geneid)
        ) +
          geom_line(alpha = 0.15, linewidth = 0.3, color = "grey50") +
          stat_summary(
            aes(group = 1),
            fun = mean,
            geom = "line",
            linewidth = 1.2,
            color = cmelo_inoculated
          ) +
          facet_wrap(~ cluster, ncol = 2) +
          labs(
            title = "Temporal clustering of inoculation-responsive genes",
            x = "Time point",
            y = "Scaled vst expression"
          ) +
          theme_classic(base_size = 12) +
          theme(
            strip.text = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.margin = margin(10, 10, 10, 10)
          )
        
        ggsave(
          file.path(figures_dir, "Cmelo_kmeans_temporal_clusters3.png"),
          cluster_plot,
          width = 12,
          height = 10,
          dpi = 300
        )
      }
    }
  }
}

############################################################
### 24. Functional enrichment analysis
###     GO / KEGG enrichment per temporal cluster
############################################################

if (!file.exists(annotation_file)) {
  
  message("Functional annotation file not found. Skipping GO/KEGG enrichment.")
  
} else if (!exists("cluster_df")) {
  
  message("Cluster assignments not available. Skipping GO/KEGG enrichment.")
  
} else if (!has_clusterprofiler || !has_enrichplot) {
  
  message("clusterProfiler and/or enrichplot are not installed. Skipping GO/KEGG enrichment.")
  
} else {
  
  annotation <- read.table(
    annotation_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  
  if (!"Geneid" %in% colnames(annotation)) {
    
    message("Annotation file lacks 'Geneid' column. Skipping enrichment.")
    
  } else {
    
    background_genes <- rownames(dds)
    cluster_ids <- sort(unique(cluster_df$cluster))
    
    ##########################################################
    ### GO enrichment
    ##########################################################
    
    if ("GO" %in% colnames(annotation)) {
      
      go_annot <- annotation %>%
        dplyr::select(Geneid, GO) %>%
        dplyr::filter(!is.na(GO), GO != "") %>%
        tidyr::separate_rows(GO, sep = ";") %>%
        dplyr::mutate(GO = trimws(GO)) %>%
        dplyr::filter(GO != "")
      
      term2gene_go <- go_annot %>%
        dplyr::select(GO, Geneid)
      
      if (nrow(term2gene_go) > 0) {
        
        for (cl in cluster_ids) {
          
          genes_cl <- cluster_df$Geneid[cluster_df$cluster == cl]
          
          run_term_enrichment(
            gene_ids = genes_cl,
            universe_ids = background_genes,
            term2gene = term2gene_go,
            out_csv = file.path(enrichment_dir, paste0("Cmelo_GO_enrichment_cluster_", cl, ".csv")),
            out_png = file.path(figures_dir, paste0("Cmelo_GO_enrichment_cluster_", cl, "_dotplot.png")),
            plot_title = paste("GO enrichment - Cluster", cl)
          )
        }
        
      } else {
        message("No valid GO mappings found. Skipping GO enrichment.")
      }
      
    } else {
      message("Annotation file lacks 'GO' column. Skipping GO enrichment.")
    }
    
    ##########################################################
    ### KEGG enrichment
    ##########################################################
    
    if ("KEGG" %in% colnames(annotation)) {
      
      kegg_annot <- annotation %>%
        dplyr::select(Geneid, KEGG) %>%
        dplyr::filter(!is.na(KEGG), KEGG != "") %>%
        tidyr::separate_rows(KEGG, sep = ";") %>%
        dplyr::mutate(KEGG = trimws(KEGG)) %>%
        dplyr::filter(KEGG != "")
      
      term2gene_kegg <- kegg_annot %>%
        dplyr::select(KEGG, Geneid)
      
      if (nrow(term2gene_kegg) > 0) {
        
        for (cl in cluster_ids) {
          
          genes_cl <- cluster_df$Geneid[cluster_df$cluster == cl]
          
          run_term_enrichment(
            gene_ids = genes_cl,
            universe_ids = background_genes,
            term2gene = term2gene_kegg,
            out_csv = file.path(enrichment_dir, paste0("Cmelo_KEGG_enrichment_cluster_", cl, ".csv")),
            out_png = file.path(figures_dir, paste0("Cmelo_KEGG_enrichment_cluster_", cl, "_dotplot.png")),
            plot_title = paste("KEGG enrichment - Cluster", cl)
          )
        }
        
      } else {
        message("No valid KEGG mappings found. Skipping KEGG enrichment.")
      }
      
    } else {
      message("Annotation file lacks 'KEGG' column. Skipping KEGG enrichment.")
    }
  }
}

############################################################
### 25. DEG summary visualization: lollipop / diverging plot
############################################################

deg_summary_condition_long <- bind_rows(
  deg_summary_condition %>%
    dplyr::select(comparison, up) %>%
    dplyr::mutate(direction = "Upregulated", count = up),
  
  deg_summary_condition %>%
    dplyr::select(comparison, down) %>%
    dplyr::mutate(direction = "Downregulated", count = -down)
) %>%
  dplyr::select(comparison, direction, count)

deg_summary_condition_long$time_point <- sub(
  "Inoculated_vs_Mock_",
  "",
  as.character(deg_summary_condition_long$comparison)
)

deg_summary_condition_long$time_point <- factor(
  deg_summary_condition_long$time_point,
  levels = time_levels
)

deg_lollipop <- ggplot(
  deg_summary_condition_long,
  aes(x = time_point, y = count, color = direction)
) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  geom_segment(
    aes(x = time_point, xend = time_point, y = 0, yend = count),
    linewidth = 1,
    alpha = 0.8
  ) +
  geom_point(size = 4) +
  coord_flip() +
  scale_y_continuous(
    labels = abs
  ) +
  scale_color_manual(
    values = deg_direction_colors
  ) +
  labs(
    x = "Comparison (Inoculated vs Mock)",
    y = "DEGs"
  ) +
  theme_cmelo_publication() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  )

ggsave(
  file.path(figures_dir, "Cmelo_DEG_lollipop_Inoculated_vs_Mock.png"),
  deg_lollipop,
  width = 8,
  height = 6,
  dpi = 300
)

############################################################
### 27. DEG lollipop:
###     Inoculated vs previous inoculated time point
############################################################

inoc_time_levels <- time_levels

if (length(inoc_time_levels) < 2) {
  
  message("Not enough time points for consecutive inoculated contrasts.")
  
} else {
  
  deg_summary_inoc_previous <- do.call(
    rbind,
    lapply(seq_along(inoc_time_levels)[-1], function(i) {
      
      tp_current <- inoc_time_levels[i]
      tp_previous <- inoc_time_levels[i - 1]
      
      res <- results(
        dds_group,
        contrast = c(
          "group",
          paste0("Inoculated_", tp_current),
          paste0("Inoculated_", tp_previous)
        )
      )
      
      comparison_label <- paste0(
        "Inoculated_",
        tp_current,
        "_vs_Inoculated_",
        tp_previous
      )
      
      output_file <- file.path(
        contrasts_dir,
        paste0("DEG_", comparison_label, ".csv")
      )
      
      filter_and_save_deg(
        res = res,
        filename = output_file,
        comparison_label = comparison_label,
        padj_cutoff = 0.05,
        lfc_cutoff = 1
      )
    })
  )
  
  write.csv(
    deg_summary_inoc_previous,
    file.path(
      contrasts_dir,
      "DEG_summary_Inoculated_vs_previous_time.csv"
    ),
    row.names = FALSE
  )
  
  ##########################################################
  ### Prepare lollipop plot data
  ##########################################################
  
  deg_summary_inoc_previous_long <- bind_rows(
    
    deg_summary_inoc_previous %>%
      dplyr::select(comparison, up) %>%
      dplyr::mutate(
        direction = "Upregulated",
        count = up
      ),
    
    deg_summary_inoc_previous %>%
      dplyr::select(comparison, down) %>%
      dplyr::mutate(
        direction = "Downregulated",
        count = -down
      )
    
  ) %>%
    dplyr::select(comparison, direction, count)
  
##########################################################
### Clean and order comparison labels
##########################################################

comparison_order <- paste0(
  inoc_time_levels[-1],
  "_vs_",
  inoc_time_levels[-length(inoc_time_levels)]
)

deg_summary_inoc_previous_long$comparison_label <- gsub(
  "Inoculated_",
  "",
  as.character(deg_summary_inoc_previous_long$comparison)
)

deg_summary_inoc_previous_long$comparison_label <- factor(
  deg_summary_inoc_previous_long$comparison_label,
  levels = comparison_order
)
  
  ##########################################################
  ### Lollipop plot
  ##########################################################
  
deg_lollipop_inoc_previous <- ggplot(
  deg_summary_inoc_previous_long,
  aes(
    x = comparison_label,
    y = count,
    color = direction
  )
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    color = "black"
  ) +
  
  geom_segment(
    aes(
      x = comparison_label,
      xend = comparison_label,
      y = 0,
      yend = count
    ),
    linewidth = 1,
    alpha = 0.8
  ) +
  
  geom_point(size = 4) +
  
  coord_flip() +
  
  scale_y_continuous(
    labels = abs
  ) +
  
  scale_color_manual(
    values = deg_direction_colors
  ) +
  
  labs(
    x = "Comparison (Sequential)",
    y = "DEGs"
  ) +
  
  theme_cmelo_publication() +
  
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  )
  
  ##########################################################
  ### Save figure
  ##########################################################
  
  ggsave(
    file.path(
      figures_dir,
      "Cmelo_DEG_lollipop_Inoculated_vs_previous_time.png"
    ),
    deg_lollipop_inoc_previous,
    width = 8,
    height = 6,
    dpi = 300
  )
}

############################################################
### 28. DEG summary for Mock vs previous Mock time point
############################################################

if (length(time_levels) < 2) {
  
  message("Not enough time points for consecutive Mock contrasts.")
  
} else {
  
  deg_summary_mock_previous <- do.call(
    rbind,
    lapply(seq_along(time_levels)[-1], function(i) {
      
      tp_current <- time_levels[i]
      tp_previous <- time_levels[i - 1]
      
      res <- results(
        dds_group,
        contrast = c(
          "group",
          paste0("Mock_", tp_current),
          paste0("Mock_", tp_previous)
        )
      )
      
      comparison_label <- paste0(
        "Mock_",
        tp_current,
        "_vs_Mock_",
        tp_previous
      )
      
      output_file <- file.path(
        contrasts_dir,
        paste0("DEG_", comparison_label, ".csv")
      )
      
      filter_and_save_deg(
        res = res,
        filename = output_file,
        comparison_label = comparison_label,
        padj_cutoff = 0.05,
        lfc_cutoff = 1
      )
    })
  )
  
  write.csv(
    deg_summary_mock_previous,
    file.path(
      contrasts_dir,
      "DEG_summary_Mock_vs_previous_time.csv"
    ),
    row.names = FALSE
  )
}

############################################################
### 26. Final message
############################################################

message("Cucumis melo analysis completed successfully.")
