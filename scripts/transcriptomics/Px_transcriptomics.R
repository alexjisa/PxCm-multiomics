############################################################
### Podosphaera xanthii RNA-seq time-course analysis
### Publication-oriented script version
###
### Includes:
### - TPM export
### - Sample distances, dendrogram, PCA, t-SNE, UMAP
### - DESeq2 across consecutive time points
### - DESeq2 across infection phases
### - Temporal clustering of dynamic genes (k-means)
### - Functional enrichment analysis (KEGG)
### - Distribution of expression changes (ordered violin plot)
###
### Main methodological decisions:
### - DESeq2 time-course analysis uses design = ~ time_point
### - Infection-phase analysis is performed separately with design = ~ infection_phase
### - Temporal clustering is performed on DEGs from consecutive comparisons
### - Expression values are averaged by time point and scaled by gene
### - k = 4 is used for k-means because it gave a better balance
###   between statistical support and biological interpretability
############################################################

############################################################
### 1. Load required libraries
############################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(pheatmap)
  library(viridis)
  library(ggrepel)
  library(Rtsne)
  library(uwot)
  library(dendextend)
  library(tidyr)
  library(stringr)
  library(grid)
  library(httr)
  library(xml2)
  library(rvest)
  library(readxl)
})

### These packages are optional and only required for enrichment analysis.
has_clusterprofiler <- requireNamespace("clusterProfiler", quietly = TRUE)
has_enrichplot <- requireNamespace("enrichplot", quietly = TRUE)
has_keggrest <- requireNamespace("KEGGREST", quietly = TRUE)

############################################################
### 2. Define input files and output directories
############################################################

counts_file <- "01_Transcriptomics_Px_counts.txt"
sampleinfo_file <- "02_Transcriptomics_Px_sample_info.txt"
annotation_file <- "Px_functional_annotation.xlsx"   # Optional

results_dir <- "results"
figures_dir <- file.path(results_dir, "figures")
tables_dir <- file.path(results_dir, "tables")
deseq_dir <- file.path(tables_dir, "deseq2")
tsne_dir <- file.path(tables_dir, "tsne")
umap_dir <- file.path(tables_dir, "umap")
clustering_dir <- file.path(tables_dir, "clustering")
enrichment_dir <- file.path(tables_dir, "enrichment")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deseq_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tsne_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(umap_dir, recursive = TRUE, showWarnings = FALSE)
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
### 4. Validate required columns
############################################################

required_sampleinfo_cols <- c("sample", "time_point", "replicate", "condition", "infection_phase")
missing_sampleinfo_cols <- setdiff(required_sampleinfo_cols, colnames(sampleinfo))

if (length(missing_sampleinfo_cols) > 0) {
  stop(
    paste(
      "The sample metadata file is missing required columns:",
      paste(missing_sampleinfo_cols, collapse = ", ")
    )
  )
}

required_counts_cols <- c("Geneid", "Length")
missing_counts_cols <- setdiff(required_counts_cols, colnames(counts))

if (length(missing_counts_cols) > 0) {
  stop(
    paste(
      "The count matrix file is missing required columns:",
      paste(missing_counts_cols, collapse = ", ")
    )
  )
}

############################################################
### 5. Validate sample names and reorder metadata
############################################################

### Assumes the first 6 columns of the count table are annotation columns
### and sample columns begin at column 7.
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

### Time points are explicitly converted into ordered factors.
### This is critical because several plots were previously being drawn
### with the time axis in the wrong order.
time_numeric <- as.integer(sub("h$", "", sampleinfo$time_point))
time_levels <- paste0(sort(unique(time_numeric)), "h")

sampleinfo$time_point <- factor(sampleinfo$time_point, levels = time_levels)
sampleinfo$replicate <- factor(sampleinfo$replicate)
sampleinfo$condition <- factor(sampleinfo$condition)
sampleinfo$infection_phase <- factor(
  sampleinfo$infection_phase,
  levels = c("Early", "Intermediate", "Late")
)

############################################################
### 7. Build count matrix and validate ordering
############################################################

count_data <- as.matrix(counts[, count_sample_names])
rownames(count_data) <- counts$Geneid
storage.mode(count_data) <- "integer"

if (!identical(colnames(count_data), rownames(sampleinfo))) {
  stop("Sample order is not identical between the count matrix and sample metadata.")
}

############################################################
### 8. Extract gene lengths and validate TPM inputs
############################################################

gene_lengths <- counts$Length
names(gene_lengths) <- counts$Geneid

if (any(is.na(gene_lengths)) || any(gene_lengths <= 0)) {
  stop("Invalid gene lengths were detected for TPM calculation.")
}


############################################################
### 9. Define annotations, colors, and reusable plot theme
############################################################

time_colors <- setNames(
  viridis::viridis(length(time_levels), option = "D"),
  time_levels
)

infection_phase_colors <- c(
  Early = viridis::viridis(12, option = "D")[12],          # purple
  Intermediate = viridis::viridis(12, option = "D")[6],   # blue-green
  Late = viridis::viridis(12, option = "D")[1]           # yellow
)

annotation_colors <- list(
  time_point = time_colors,
  infection_phase = infection_phase_colors
)

sample_annotation <- sampleinfo[, c("time_point", "infection_phase"), drop = FALSE]

heatmap_colors <- viridis::viridis(255, option = "D")

theme_px_publication <- function(base_size = 12) {
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

### TPM is exported because it is useful for downstream exploration,
### visualization and biological interpretation, although DESeq2 uses raw counts.
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
  file.path(tables_dir, "Px_counts_TPM.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

############################################################
### 11. Build DESeq2 dataset for time-point analysis
############################################################

### Time-course DESeq2 model:
### design = ~ time_point
###
### Justification:
### time_point provides the finest temporal resolution.
### infection_phase is analyzed separately because it is biologically derived
### from time and should not be modeled together here if it is strongly confounded.
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sampleinfo,
  design = ~ time_point
)

############################################################
### 12. Filter low-expression genes and normalize counts
############################################################

### Genes are kept if they have at least 10 counts in at least 3 samples.
### This reduces noise before transformation and DE analysis.
keep_genes <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep_genes, ]

dds <- estimateSizeFactors(dds)

message(sum(keep_genes), " genes retained after low-count filtering.")

############################################################
### 13. Apply rlog transformation for exploratory analyses
############################################################

### rlog is used here for exploratory analyses because it stabilizes variance
### and is appropriate for PCA, clustering and distance-based plots.
rld <- rlog(dds, blind = FALSE)

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

filter_and_save_deg <- function(res, filename, comparison_label, padj_cutoff = 0.05, lfc_cutoff = 1) {
  res_df <- as.data.frame(res)
  res_df$Geneid <- rownames(res_df)
  res_df <- res_df[, c("Geneid", setdiff(colnames(res_df), "Geneid"))]
  
  sig_deg <- res_df[
    !is.na(res_df$padj) &
      res_df$padj < padj_cutoff &
      abs(res_df$log2FoldChange) > lfc_cutoff,
  ]
  
  sig_deg_up <- sig_deg[
    sig_deg$log2FoldChange > lfc_cutoff,
  ]
  
  sig_deg_down <- sig_deg[
    sig_deg$log2FoldChange < -lfc_cutoff,
  ]
  
  up <- nrow(sig_deg_up)
  down <- nrow(sig_deg_down)
  
  write.csv(sig_deg, filename, row.names = FALSE)
  
  base_name <- tools::file_path_sans_ext(basename(filename))
  output_dir <- dirname(filename)
  
  write.table(
    sig_deg_up$Geneid,
    file = file.path(output_dir, paste0(base_name, "_UP_genes.txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  write.table(
    sig_deg_down$Geneid,
    file = file.path(output_dir, paste0(base_name, "_DOWN_genes.txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  data.frame(
    comparison = comparison_label,
    total = nrow(sig_deg),
    up = up,
    down = down,
    stringsAsFactors = FALSE
  )
}

run_term_enrichment <- function(
    gene_ids,
    universe_ids,
    term2gene,
    term2name = NULL,
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
  
  if (is.null(term2name)) {
    enr <- clusterProfiler::enricher(
      gene = gene_ids,
      universe = universe_ids,
      TERM2GENE = term2gene,
      minGSSize = 1,
      maxGSSize = 10000,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  } else {
    enr <- clusterProfiler::enricher(
      gene = gene_ids,
      universe = universe_ids,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      minGSSize = 1,
      maxGSSize = 10000,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05
    )
  }
  
  if (is.null(enr)) {
    return(NULL)
  }
  
  enr_df <- as.data.frame(enr)
  if (nrow(enr_df) == 0) {
    return(NULL)
  }

  if (!is.null(term2name) && nrow(term2name) > 0 && all(c("ID", "Description") %in% names(enr_df))) {
    term_lookup <- term2name %>%
      dplyr::distinct(term, .keep_all = TRUE)

    resolved_from_id <- term_lookup$name[match(enr_df$ID, term_lookup$term)]
    resolved_from_desc <- term_lookup$name[match(enr_df$Description, term_lookup$term)]

    enr_df$Description <- dplyr::coalesce(
      resolved_from_id,
      resolved_from_desc,
      enr_df$Description
    )
  }

  unresolved_kegg_terms <- grepl("^(ko|map)\\d{5}$", enr_df$Description)
  if (any(unresolved_kegg_terms)) {
    message(
      "Dropping unresolved KEGG terms without names: ",
      paste(unique(enr_df$Description[unresolved_kegg_terms]), collapse = ", ")
    )
    enr_df <- enr_df[!unresolved_kegg_terms, , drop = FALSE]
  }

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
### 15. Generate sample-to-sample distance heatmap & dendrogram
############################################################

sample_dists_euclidean <- dist(t(assay(rld)), method = "euclidean")
sample_dist_matrix_euclidean <- as.matrix(sample_dists_euclidean)
rownames(sample_dist_matrix_euclidean) <- colnames(rld)
colnames(sample_dist_matrix_euclidean) <- colnames(rld)

save_heatmap(
  sample_dist_matrix_euclidean,
  sample_dists_euclidean,
  heatmap_colors,
  file.path(figures_dir, "Px_sample_distance_heatmap_euclidean.png")
)

sample_dists <- dist(t(assay(rld)), method = "euclidean")
hc <- hclust(sample_dists, method = "ward.D2")
dend <- as.dendrogram(hc)

labels_colors(dend) <- infection_phase_colors[
  as.character(sampleinfo$infection_phase)
][order.dendrogram(dend)]

dend <- set(dend, "labels_cex", 0.9)
max_h <- attr(dend, "height")

png(
  file.path(figures_dir, "Px_sample_dendrogram_infection_phase.png"),
  width = 2400,
  height = 2200,
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
### 16. Perform PCA and export clean plot
############################################################

pca_data <- plotPCA(
  rld,
  intgroup = "infection_phase",
  returnData = TRUE
)

pca_data$infection_phase <- factor(
  pca_data$infection_phase,
  levels = c("Early", "Intermediate", "Late")
)

pca_data$sample <- rownames(pca_data)

percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = infection_phase,
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
  scale_color_manual(
    name = "Infection phase",
    values = infection_phase_colors
  ) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_y_continuous(expand = expansion(mult = 0.02)) +
  labs(
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  ) +
  theme_px_publication()

ggsave(
  file.path(figures_dir, "Px_PCA_infection_phase.png"),
  pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

############################################################
### 17. Run t-SNE with all samples
############################################################

set.seed(42)

n_samples <- ncol(rld)
tsne_perplexity <- min(3, floor((n_samples - 1) / 3))

if (tsne_perplexity < 1) {
  
  message("Not enough samples for t-SNE. Skipping t-SNE.")
  
} else {
  
  tsne_out <- Rtsne(
    t(assay(rld)),
    perplexity = tsne_perplexity,
    check_duplicates = FALSE
  )
  
  tsne_data <- data.frame(
    tSNE1 = tsne_out$Y[, 1],
    tSNE2 = tsne_out$Y[, 2],
    sample = colnames(rld),
    infection_phase = factor(
      sampleinfo$infection_phase,
      levels = c("Early", "Intermediate", "Late")
    )
  )
  
  write.csv(
    tsne_data,
    file.path(tsne_dir, "Px_tSNE_coordinates.csv"),
    row.names = FALSE
  )
  
  tsne_plot <- ggplot(
    tsne_data,
    aes(
      x = tSNE1,
      y = tSNE2,
      color = infection_phase,
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
    scale_color_manual(
      name = "Infection phase",
      values = infection_phase_colors
    ) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(
      x = "t-SNE 1",
      y = "t-SNE 2"
    ) +
    theme_px_publication()
  
  ggsave(
    file.path(figures_dir, "Px_tSNE_infection_phase_labeled.png"),
    tsne_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

############################################################
### 18. Run UMAP with all samples
############################################################

set.seed(42)

if (n_samples < 3) {
  
  message("Not enough samples for UMAP. Skipping UMAP.")
  
} else {
  
  umap_n_neighbors <- min(15, n_samples - 1)
  umap_n_neighbors <- max(2, umap_n_neighbors)
  
  umap_out <- uwot::umap(
    t(assay(rld)),
    n_neighbors = umap_n_neighbors,
    min_dist = 0.3,
    metric = "euclidean",
    verbose = FALSE,
    ret_model = FALSE
  )
  
  umap_data <- data.frame(
    UMAP1 = umap_out[, 1],
    UMAP2 = umap_out[, 2],
    sample = colnames(rld),
    infection_phase = factor(
      sampleinfo$infection_phase,
      levels = c("Early", "Intermediate", "Late")
    )
  )
  
  write.csv(
    umap_data,
    file.path(umap_dir, "Px_UMAP_coordinates.csv"),
    row.names = FALSE
  )
  
  umap_plot <- ggplot(
    umap_data,
    aes(
      x = UMAP1,
      y = UMAP2,
      color = infection_phase,
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
    scale_color_manual(
      name = "Infection phase",
      values = infection_phase_colors
    ) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_px_publication()
  
  ggsave(
    file.path(figures_dir, "Px_UMAP_infection_phase_labeled.png"),
    umap_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

############################################################
### 19. Run DESeq2 differential expression analysis (time points)
############################################################

dds <- DESeq(dds)

############################################################
### 20. Define consecutive time-point comparisons
############################################################

### Consecutive comparisons preserve temporal continuity and allow
### the identification of transcriptomic shifts between adjacent stages.
time_levels_ordered <- levels(droplevels(sampleinfo$time_point))

comparison_pairs <- data.frame(
  current_time = time_levels_ordered[-1],
  previous_time = time_levels_ordered[-length(time_levels_ordered)],
  stringsAsFactors = FALSE
)

############################################################
### 21. Export DEGs for each time point versus the previous one
############################################################

deg_summary <- do.call(
  rbind,
  lapply(seq_len(nrow(comparison_pairs)), function(i) {
    current_time <- comparison_pairs$current_time[i]
    previous_time <- comparison_pairs$previous_time[i]
    
    res <- results(dds, contrast = c("time_point", current_time, previous_time))
    
    comparison_label <- paste0(current_time, "_vs_", previous_time)
    output_file <- file.path(deseq_dir, paste0("DEG_", comparison_label, ".csv"))
    
    filter_and_save_deg(res, output_file, comparison_label)
  })
)

write.csv(
  deg_summary,
  file.path(deseq_dir, "DEG_summary_consecutive_timepoints.csv"),
  row.names = FALSE
)

############################################################
### 22. Plot number of up- and downregulated genes
###     between consecutive time points
############################################################

deg_summary_long <- bind_rows(
  deg_summary %>%
    dplyr::select(comparison, up) %>%
    dplyr::mutate(direction = "Upregulated", count = up),
  deg_summary %>%
    dplyr::select(comparison, down) %>%
    dplyr::mutate(direction = "Downregulated", count = -down)
) %>%
  dplyr::select(comparison, direction, count)

### Force chronological order
comparison_order <- paste0(
  comparison_pairs$current_time,
  "_vs_",
  comparison_pairs$previous_time
)

deg_summary_long$comparison <- factor(
  deg_summary_long$comparison,
  levels = comparison_order
)

############################################################
### 22.1 Original barplot
############################################################

deg_barplot <- ggplot(
  deg_summary_long,
  aes(x = comparison, y = count, fill = direction)
) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(
    values = c(
      Upregulated = "#4575B4",
      Downregulated = "#D73027"
    )
  ) +
  labs(
    title = "Differentially expressed genes between consecutive time points",
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
  file.path(figures_dir, "Px_DEG_summary_consecutive_timepoints.png"),
  deg_barplot,
  width = 9,
  height = 7,
  dpi = 300
)

############################################################
### 21.2 Lollipop / diverging plot
############################################################

deg_summary_long$comparison_clean <- as.character(deg_summary_long$comparison)

deg_summary_long$comparison_clean <- factor(
  deg_summary_long$comparison_clean,
  levels = comparison_order
)

deg_lollipop <- ggplot(
  deg_summary_long,
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
  scale_y_continuous(
    labels = abs
  ) +
  scale_color_manual(
    values = c(
      Upregulated = "#D62728",
      Downregulated = "#1F77B4"
    )
  ) +
  labs(
    x = "Comparison (Sequential)",
    y = "DEGs"
  ) +
  theme_px_publication() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  )

ggsave(
  file.path(figures_dir, "Px_DEG_lollipop_consecutive_timepoints.png"),
  deg_lollipop,
  width = 8,
  height = 6,
  dpi = 300
)
############################################################
### 22. Differential expression across infection phases
############################################################

### Phase-based DE analysis is run separately because infection_phase is
### a biologically meaningful categorical summary of the trajectory.
dds_phase <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sampleinfo,
  design = ~ infection_phase
)

keep_genes_phase <- rowSums(counts(dds_phase) >= 10) >= 3
dds_phase <- dds_phase[keep_genes_phase, ]
dds_phase <- DESeq(dds_phase)

phase_comparisons <- list(
  c("infection_phase", "Intermediate", "Early"),
  c("infection_phase", "Late", "Intermediate"),
  c("infection_phase", "Late", "Early")
)

deg_summary_phase <- do.call(
  rbind,
  lapply(phase_comparisons, function(comp) {
    res <- results(dds_phase, contrast = comp)
    comparison_label <- paste0(comp[2], "_vs_", comp[3])
    output_file <- file.path(deseq_dir, paste0("DEG_phase_", comparison_label, ".csv"))
    filter_and_save_deg(res, output_file, comparison_label)
  })
)

write.csv(
  deg_summary_phase,
  file.path(deseq_dir, "DEG_summary_infection_phases.csv"),
  row.names = FALSE
)

deg_summary_phase_long <- bind_rows(
  deg_summary_phase %>%
    dplyr::select(comparison, up) %>%
    dplyr::mutate(direction = "Upregulated", count = up),
  deg_summary_phase %>%
    dplyr::select(comparison, down) %>%
    dplyr::mutate(direction = "Downregulated", count = -down)
) %>%
  dplyr::select(comparison, direction, count)

phase_barplot <- ggplot(deg_summary_phase_long, aes(x = comparison, y = count, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c(Upregulated = "#4575B4", Downregulated = "#D73027")) +
  labs(
    title = "Differentially expressed genes across infection phases",
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
  file.path(figures_dir, "Px_DEG_summary_infection_phases.png"),
  phase_barplot,
  width = 8,
  height = 6,
  dpi = 300
)

############################################################
### 23. Distribution of expression changes (ordered violin plot)
############################################################

### This figure is kept as a complementary summary of global LFC distributions.
### It is made more informative by:
### - enforcing chronological order
### - adding the number of significant DEGs on top of each violin
### - restricting the visible y-range so extreme outliers do not flatten the plot
if (!exists("comparison_pairs") || nrow(comparison_pairs) == 0) {
  message("comparison_pairs not found or empty. Skipping violin plot.")
} else {
  
  lfc_list <- lapply(seq_len(nrow(comparison_pairs)), function(i) {
    current_time <- comparison_pairs$current_time[i]
    previous_time <- comparison_pairs$previous_time[i]
    
    res <- results(dds, contrast = c("time_point", current_time, previous_time))
    res_df <- as.data.frame(res)
    res_df$Geneid <- rownames(res_df)
    res_df$comparison <- paste0(current_time, "_vs_", previous_time)
    res_df
  })
  
  lfc_data <- bind_rows(lfc_list)
  
  if (nrow(lfc_data) == 0) {
    message("No differential expression results available for violin plot.")
  } else {
    
    lfc_data <- lfc_data %>%
      dplyr::filter(!is.na(log2FoldChange))
    
    if (nrow(lfc_data) == 0) {
      message("No valid log2FoldChange values available for violin plot.")
    } else {
      
      comparison_order <- paste0(
        comparison_pairs$current_time,
        "_vs_",
        comparison_pairs$previous_time
      )
      
      lfc_data$comparison <- factor(
        lfc_data$comparison,
        levels = comparison_order
      )
      
      write.csv(
        lfc_data,
        file.path(deseq_dir, "LFC_all_consecutive_timepoints.csv"),
        row.names = FALSE
      )
      
      deg_counts_labels <- lfc_data %>%
        dplyr::mutate(
          is_deg = !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1
        ) %>%
        dplyr::group_by(comparison) %>%
        dplyr::summarise(
          n_deg = sum(is_deg, na.rm = TRUE),
          .groups = "drop"
        )
      
      deg_counts_labels$label <- paste0("n=", deg_counts_labels$n_deg)
      
      lfc_violin <- ggplot(lfc_data, aes(x = comparison, y = log2FoldChange)) +
        geom_violin(fill = "#BFDCEB", color = "black", trim = TRUE, linewidth = 0.4) +
        geom_boxplot(
          width = 0.12,
          outlier.shape = NA,
          fill = "white",
          color = "black",
          linewidth = 0.4
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.6) +
        geom_text(
          data = deg_counts_labels,
          aes(x = comparison, y = 3.6, label = label),
          inherit.aes = FALSE,
          size = 3.2,
          fontface = "bold"
        ) +
        coord_cartesian(ylim = c(-4, 4)) +
        labs(
          title = "Distribution of expression changes across consecutive time points",
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
      
      print(lfc_violin)
      
      ggsave(
        file.path(figures_dir, "Px_log2FC_violin_consecutive_timepoints.png"),
        plot = lfc_violin,
        width = 11,
        height = 6.5,
        dpi = 300
      )
      
      message("Violin plot saved successfully.")
    }
  }
}

############################################################
### 24. Temporal clustering of dynamic genes (k-means)
############################################################

### Genes for temporal clustering are defined as the union of DEGs detected
### in consecutive time-point comparisons.
###
### Justification:
### this focuses the clustering on genes that actually change during the time course,
### instead of diluting patterns with mostly static genes.
deg_files_time <- list.files(
  deseq_dir,
  pattern = "^DEG_[0-9]+h_vs_[0-9]+h\\.csv$",
  full.names = TRUE
)

if (length(deg_files_time) == 0) {
  message("No DEG files found for consecutive time-point comparisons. Skipping clustering.")
} else {
  
  deg_lists_time <- lapply(deg_files_time, read.csv)
  deg_genes_time <- unique(unlist(lapply(deg_lists_time, function(x) x$Geneid)))
  deg_genes_time <- intersect(deg_genes_time, rownames(rld))
  
  if (length(deg_genes_time) < 2) {
    message("Not enough dynamic genes for k-means clustering. Skipping clustering.")
  } else {
    
    rld_mat <- assay(rld)
    rld_deg <- rld_mat[deg_genes_time, , drop = FALSE]
    
    ### Replicates are averaged by time point to obtain one temporal profile per gene.
    ### This emphasizes trajectory shape rather than replicate-level variation.
    time_means <- sapply(levels(sampleinfo$time_point), function(tp) {
      rowMeans(rld_deg[, sampleinfo$time_point == tp, drop = FALSE])
    })
    
    time_means <- as.matrix(time_means)
    colnames(time_means) <- levels(sampleinfo$time_point)
    
    ### Values are scaled by gene so clustering captures temporal pattern,
    ### not absolute expression magnitude.
    time_means_scaled <- t(scale(t(time_means)))
    time_means_scaled <- time_means_scaled[complete.cases(time_means_scaled), , drop = FALSE]
    
    n_genes_cluster <- nrow(time_means_scaled)
    
    ### k = 4 was selected after comparing k values using elbow and silhouette criteria.
    ### It produced a cleaner and more interpretable solution than k = 5 or k = 6,
    ### while still preserving the main biological patterns.
    k <- min(4, n_genes_cluster)
    
    if (k < 2) {
      message("Not enough genes remaining after scaling for k-means clustering. Skipping.")
    } else {
      
      set.seed(123)
      km <- kmeans(time_means_scaled, centers = k, nstart = 50)
      
      cluster_df <- as.data.frame(time_means_scaled)
      cluster_df$Geneid <- rownames(time_means_scaled)
      cluster_df$cluster <- factor(km$cluster)
      
      write.csv(
        cluster_df[, c("Geneid", "cluster")],
        file.path(clustering_dir, "Px_temporal_gene_clusters.csv"),
        row.names = FALSE
      )
      
      cluster_centers <- as.data.frame(km$centers)
      cluster_centers$cluster <- factor(seq_len(nrow(cluster_centers)))
      
      write.csv(
        cluster_centers,
        file.path(clustering_dir, "Px_temporal_cluster_centers.csv"),
        row.names = FALSE
      )
      
      cluster_long <- cluster_df %>%
        tidyr::pivot_longer(
          cols = all_of(levels(sampleinfo$time_point)),
          names_to = "time_point",
          values_to = "scaled_expression"
        )
      
      cluster_long$time_point <- factor(
        cluster_long$time_point,
        levels = levels(sampleinfo$time_point)
      )
      
      ### The facet layout is explicitly set to 2 columns so the 4-cluster solution
      ### is exported as a clean 2x2 panel figure without truncation.
      cluster_plot <- ggplot(cluster_long, aes(x = time_point, y = scaled_expression, group = Geneid)) +
        geom_line(alpha = 0.15, linewidth = 0.3, color = "grey50") +
        stat_summary(
          aes(group = 1),
          fun = mean,
          geom = "line",
          linewidth = 1.2,
          color = "red"
        ) +
        facet_wrap(~ cluster, ncol = 2) +
        labs(
          title = "Temporal clustering of differentially expressed genes",
          x = "Time point",
          y = "Scaled rlog expression"
        ) +
        theme_classic(base_size = 12) +
        theme(
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(10, 10, 10, 10)
        )
      
      ggsave(
        file.path(figures_dir, "Px_kmeans_temporal_clusters.png"),
        cluster_plot,
        width = 12,
        height = 10,
        dpi = 300
      )
    }
  }
}

############################################################
### 25. Functional enrichment analysis
###     KEGG enrichment per consecutive time-point comparison
############################################################

kegg_dir <- file.path(enrichment_dir, "KEGG_eggNOG_results")
kegg_fig_dir <- file.path(kegg_dir, "figures")

dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(kegg_fig_dir, recursive = TRUE, showWarnings = FALSE)

clean_gene <- function(x) {
  x <- str_trim(x)
  x[!is.na(x) & x != "" & x != "-"] %>% unique()
}

if (!file.exists(annotation_file)) {
  
  message("Functional annotation file not found. Skipping KEGG enrichment.")
  
} else if (!exists("comparison_pairs") || nrow(comparison_pairs) == 0) {
  
  message("Comparison pairs not available. Skipping KEGG enrichment.")
  
} else {
  
  annotation <- readxl::read_excel(annotation_file) %>%
    as.data.frame()
  
  if (!"Geneid" %in% colnames(annotation)) {
    message("Annotation file lacks 'Geneid' column. Skipping KEGG enrichment.")
  } else {
    
    comparison_order <- paste0(
      comparison_pairs$current_time,
      "_vs_",
      comparison_pairs$previous_time
    )

    if (!"KEGG_ko" %in% colnames(annotation)) {
      message("Annotation file lacks 'KEGG_ko' column. Skipping KEGG enrichment.")
    } else {
      kegg_annot <- annotation %>%
        dplyr::select(Geneid, KEGG_ko) %>%
        dplyr::filter(
          !is.na(KEGG_ko),
          KEGG_ko != "-",
          KEGG_ko != ""
        ) %>%
        tidyr::separate_rows(KEGG_ko, sep = ",") %>%
        dplyr::mutate(KEGG_ko = trimws(KEGG_ko)) %>%
        dplyr::filter(KEGG_ko != "") %>%
        dplyr::mutate(KEGG_ko = gsub("^ko:", "", KEGG_ko))

      term2gene_kegg <- kegg_annot %>%
        dplyr::select(KEGG_ko, Geneid)

      if (nrow(term2gene_kegg) == 0) {
        message("No valid KEGG mappings found. Skipping KEGG enrichment.")
      } else {

        term2name_kegg <- NULL
        if (has_keggrest) {
          ko_list <- tryCatch(KEGGREST::keggList("ko"), error = function(e) character())
          if (length(ko_list) > 0) {
            term2name_kegg <- tibble::tibble(
              term = sub("^ko:", "", names(ko_list)),
              name = unname(ko_list)
            ) %>%
              dplyr::filter(term %in% term2gene_kegg$KEGG_ko)
          }
        }

        summary_list <- list()

        for (i in seq_len(nrow(comparison_pairs))) {
          current_time <- comparison_pairs$current_time[i]
          previous_time <- comparison_pairs$previous_time[i]
          comparison_label <- paste0(current_time, "_vs_", previous_time)

          for (direction in c("UP", "DOWN")) {
            genes_file <- file.path(
              deseq_dir,
              paste0("DEG_", comparison_label, "_", direction, "_genes.txt")
            )

            if (!file.exists(genes_file)) {
              next
            }

            gene_ids <- readLines(genes_file, warn = FALSE) %>%
              clean_gene()

            if (length(gene_ids) < 2) {
              next
            }

            message(
              "Procesando KEGG para ",
              comparison_label,
              " ",
              direction,
              " (",
              length(gene_ids),
              " genes)"
            )

            enr_df <- tryCatch(
              run_term_enrichment(
                gene_ids = gene_ids,
                universe_ids = rownames(dds),
                term2gene = term2gene_kegg,
                term2name = term2name_kegg,
                out_csv = file.path(
                  kegg_dir,
                  paste0("DEG_", comparison_label, "_", direction, "_KEGG_ko.csv")
                ),
                out_png = file.path(
                  kegg_fig_dir,
                  paste0("DEG_", comparison_label, "_", direction, "_KEGG_ko.pdf")
                ),
                plot_title = paste("KEGG KO enrichment -", comparison_label, direction)
              ),
              error = function(e) {
                message("Error KEGG para ", genes_file, ": ", e$message)
                return(NULL)
              }
            )

            if (is.null(enr_df) || nrow(enr_df) == 0) {
              next
            }

            filtered <- enr_df %>%
              mutate(
                pathway = Description
              ) %>%
              mutate(
                minus_log10_FDR = -log10(p.adjust),
                hour = as.numeric(str_match(comparison_label, "^(\\d+)h")[, 2]),
                direction = direction,
                comparison = paste0("DEG_", comparison_label, "_", direction),
                comparison_label = comparison_label
              ) %>%
              dplyr::select(
                pathway,
                minus_log10_FDR,
                hour,
                direction,
                comparison,
                comparison_label,
                Count,
                pvalue,
                p.adjust,
                qvalue,
                geneID
              )

            summary_list[[paste0("DEG_", comparison_label, "_", direction)]] <- filtered
          }
        }

        if (length(summary_list) > 0) {
          all_results <- bind_rows(summary_list) %>%
            filter(!grepl("^(ko|map)\\d{5}$", pathway)) %>%
            arrange(hour, direction, p.adjust)

          all_results <- all_results %>%
            mutate(
              comparison_label = factor(comparison_label, levels = comparison_order),
              direction = factor(direction, levels = c("DOWN", "UP"))
            )

          write_tsv(
            all_results,
            file.path(kegg_dir, "ALL_KEGG_results.tsv")
          )

          excluded_pathways <- c(
            "Metabolic pathways",
            "Biosynthesis of secondary metabolites"
          )

          top_pathways <- all_results %>%
            filter(!pathway %in% excluded_pathways) %>%
            group_by(pathway) %>%
            summarise(
              n_sig = n(),
              best_score = max(minus_log10_FDR, na.rm = TRUE),
              mean_score = mean(minus_log10_FDR, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            arrange(desc(n_sig), desc(best_score), desc(mean_score), pathway) %>%
            slice_head(n = 16)

          plot_data <- all_results %>%
            semi_join(top_pathways, by = "pathway") %>%
            mutate(
              comparison_label = factor(comparison_label, levels = comparison_order),
              direction = factor(direction, levels = c("DOWN", "UP"))
            )

          if (nrow(plot_data) > 0) {
            pathway_levels <- top_pathways %>%
              arrange(desc(n_sig), desc(best_score), desc(mean_score), pathway) %>%
              pull(pathway)

            plot_data <- plot_data %>%
              mutate(
                pathway_plot = factor(pathway, levels = pathway_levels)
              )

            p_global <- ggplot(
              plot_data,
              aes(
                x = comparison_label,
                y = pathway_plot,
                fill = minus_log10_FDR
              )
            ) +
              geom_tile() +
              scale_fill_gradient(low = "grey98", high = "grey60") +
              facet_grid(. ~ direction) +
              scale_x_discrete(drop = FALSE) +
              scale_y_discrete(drop = FALSE) +
              labs(
                x = NULL,
                y = NULL,
                fill = "-log10(FDR)"
              ) +
              theme_minimal(base_size = 11) +
              theme(
                text = element_text(size = 11),
                legend.position = "right",
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 9),
                strip.text = element_text(size = 10, face = "bold"),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
                axis.text.y = element_text(size = 9, face = "bold"),
                axis.ticks = element_blank(),
                panel.grid = element_blank(),
                panel.spacing.x = grid::unit(0, "pt"),
                panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
              )

            ggsave(
              file.path(kegg_fig_dir, "ALL_KEGG_timecourse_heatmap.pdf"),
              p_global,
              width = 16,
              height = 7.5
            )
          }
        }
      }
    }
  }
}
############################################################
### 26. Final message
############################################################

message("Analysis completed successfully.")
