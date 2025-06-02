#' Find potential markers using DE and PCC analysis
#' @title findMechanismMarkers
#' @description Performs differential expression (DE) and Pearson correlation coefficient (PCC) analysis to identify potential marker genes associated with a user-defined cellular mechanism (e.g., maturation, proliferation) across two populations.
#'
#' @param seurat_obj A labeled Seurat object. The Seurat object returned from isletAutoLabeler() is ready for use.
#' @param ident.1 Character. Cell identity label for population 1.
#' @param ident.2 Character. Cell identity label for population 2.
#' @param top_n_DEGs_volcano Integer. Number of DEGs to visualize in the volcano plot (default = 20).
#' @param top_n_DEGs_pcc Integer. Number of top DE genes to consider in PCC (default = 500).
#' @param direction_DEGs Character. "up" or "down" to select direction of DE genes (default = NA, no direction filtering).
#' @param logfc.threshold Numeric. Log2FC threshold for DE filtering (default = 0).
#' @param pval.threshold Numeric. p-value threshold for DE filtering (default = 0.05).
#' @param marker_genes Character vector. Marker genes to calculate PCC against (input required).
#' @param corr_threshold Numeric. Absolute R threshold for selecting strongly correlated genes (default = 0.5).
#' @return A list containing:
#'   - full_DEGs: full DEGs from DE analysis
#'   - top_DEGs: DEGs used for volcano plot (not filtered by direction)
#'   - top_DEGs_filtered: DEGs used for PCC (filtered by direction and LogFC thresholds)
#'   - volcano_plot: ggplot object of volcano plot
#'   - cor_matrix_raw: full PCC matrix (raw, unfiltered)
#'   - cor_matrix_filtered: PCC matrix (cells expressing both genes)
#'   - strongly_correlated_genes: genes with abs(R) > threshold using raw data
#'   - strongly_correlated_genes_filtered: genes with abs(R) > threshold using filtered data
#'   - markerwise_correlation_tables_raw: correlation tables per marker (raw)
#'   - combined_table_raw: gene-wise PCC with all marker genes (raw)
#'   - markerwise_correlation_tables_filtered: correlation tables per marker
#'   - combined_table_filtered: gene-wise PCC with all marker genes
#' @import Seurat, ggplot2, pheatmap
#' @export

findMechanismMarkers <- function(
    seurat_obj,
    ident.1,
    ident.2,
    top_n_DEGs_volcano = 20,
    top_n_DEGs_pcc = 500,
    direction_DEGs = NA,
    logfc.threshold = 0,
    pval.threshold = 0.05,
    marker_genes,
    corr_threshold = 0.5
) {
  library(ggplot2)
  library(dplyr)

  Idents(seurat_obj) <- "CellType"
  seurat_obj <- subset(seurat_obj, idents = c(ident.1, ident.2))
  seurat_obj <- JoinLayers(seurat_obj)

  # DE for visualization (volcano plot)
  DEG_full <- FindMarkers(seurat_obj, ident.1 = ident.1, ident.2 = ident.2)
  volcano_genes <- head(DEG_full[order(DEG_full$p_val), ], top_n_DEGs_volcano)
  volcano_genes$diffexpressed <- ifelse(
    volcano_genes$avg_log2FC > 0, "Upregulated",
    ifelse(volcano_genes$avg_log2FC < 0, "Downregulated", "Not significant")
  )
  volcano_genes$delabel <- rownames(volcano_genes)
  volcano_plot <- ggplot(volcano_genes, aes(x = avg_log2FC, y = -log10(p_val), color = diffexpressed, label = delabel)) +
    geom_point(size = 1.2) +
    geom_text_repel(max.overlaps = 20, size = 3) +
    scale_color_manual(values = c("Downregulated" = "blue", "Not significant" = "grey", "Upregulated" = "red")) +
    labs(x = "log2 Fold Change", y = "-log10 p-value", title = paste("Volcano plot:", ident.1, "vs", ident.2))

  # Subset DEGs for PCC; optional filtering by direction
  DEG_filtered <- DEG_full[DEG_full$p_val < pval.threshold, ]
  if (!is.na(direction_DEGs)) {
    if (direction_DEGs == "down") {
      DEG_filtered <- subset(DEG_filtered, avg_log2FC < logfc.threshold)
    } else if (direction_DEGs == "up") {
      DEG_filtered <- subset(DEG_filtered, avg_log2FC > logfc.threshold)
    }
  }
  DEG_filtered <- DEG_filtered[order(DEG_filtered$avg_log2FC), ]
  top_genes <- head(rownames(DEG_filtered), top_n_DEGs_pcc)
  top_genes <- unique(c(top_genes, marker_genes))

  expr_data <- as.matrix(GetAssayData(seurat_obj, slot = "data")[top_genes, ])
  cor_matrix_raw <- cor(t(expr_data), method = "pearson")

  # Raw correlation (all cells)
  strongly_correlated_raw <- list()
  per_marker_tables_raw <- list()
  for (marker in marker_genes) {
    r_vals_raw <- cor_matrix_raw[marker, top_genes]
    sig_genes_raw <- names(r_vals_raw[!is.na(r_vals_raw) & abs(r_vals_raw) > corr_threshold])
    strongly_correlated_raw[[marker]] <- sig_genes_raw
    per_marker_tables_raw[[marker]] <- data.frame(Gene = names(r_vals_raw), PearsonR = r_vals_raw)
  }
  # Combined table of R values per marker
  valid_tables_raw <- per_marker_tables_raw[sapply(per_marker_tables_raw, function(df) is.data.frame(df) && "Gene" %in% colnames(df))]
  combined_table_raw <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), valid_tables_raw)

  # Filtered correlation (co-expressed cells only)
  cor_matrix_filtered <- matrix(NA, nrow = length(top_genes), ncol = length(top_genes), dimnames = list(top_genes, top_genes))
  for (i in seq_along(top_genes)) {
    for (j in seq_along(top_genes)) {
      gene1 <- top_genes[i]
      gene2 <- top_genes[j]
      valid_cells <- expr_data[gene1, ] > 0 & expr_data[gene2, ] > 0
      if (sum(valid_cells) > 2) {
        cor_matrix_filtered[gene1, gene2] <- cor(expr_data[gene1, valid_cells], expr_data[gene2, valid_cells])
      }
    }
  }

  strongly_correlated_filtered <- list()
  per_marker_tables <- list()
  for (marker in marker_genes) {
    r_vals_filtered <- cor_matrix_filtered[marker, top_genes]
    sig_genes <- names(r_vals_filtered[!is.na(r_vals_filtered) & abs(r_vals_filtered) > corr_threshold])
    strongly_correlated_filtered[[marker]] <- sig_genes
    per_marker_tables[[marker]] <- data.frame(Gene = names(r_vals_filtered), PearsonR = r_vals_filtered)
  }
  valid_tables_filtered <- per_marker_tables[sapply(per_marker_tables, function(df) is.data.frame(df) && "Gene" %in% colnames(df))]
  combined_table <- Reduce(function(x, y) merge(x, y, by = "Gene", all = TRUE), valid_tables_filtered)

  return(list(
    full_DEGs = DEG_full,
    top_DEGs = volcano_genes,
    top_DEGs_filtered = DEG_filtered,
    volcano_plot = volcano_plot,
    cor_matrix_raw = cor_matrix_raw,
    cor_matrix_filtered = cor_matrix_filtered,
    strongly_correlated_genes = strongly_correlated_raw,
    strongly_correlated_genes_filtered = strongly_correlated_filtered,
    markerwise_correlation_tables_raw = per_marker_tables_raw,
    combined_table_raw = combined_table_raw,
    markerwise_correlation_tables_filtered = per_marker_tables,
    combined_table_filtered = combined_table
  ))
}

