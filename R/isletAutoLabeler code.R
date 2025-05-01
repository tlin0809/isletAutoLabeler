#' Auto-label islet cell types
#'
#' Automatically labels Alpha, Beta, Bihormonal, and Other cell clusters in a Seurat object using average marker expression.
#' Example Usage: isletAutoLabeler(seurat_obj, assay = "RNA", cluster_col = "unintegrated_clusters", source = "SC")
#'
#' @param seurat_obj A Seurat object to be labeled.
#' @param assay The assay to use (default = "integrated")
#' @param cluster_col Column in metadata to use for clustering (default = "seurat_clusters")
#' @param source Cell source, used as a prefix for the cell type label (e.g., "Human", "SC")
#' @return Updated Seurat object with a new column named "CellType".
#' @export

isletAutoLabeler <- function(seurat_obj, assay = "integrated",
                             cluster_col = "seurat_clusters", source = "Human") {

  # Set prefix for cell type label
  prefix <- source

  # Extract expression matrices
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, slot = "data")

  # Calculate global average expression for each marker
  mean_GCG  <- mean(expr_matrix["GCG", ])
  mean_ARX  <- mean(expr_matrix["ARX", ])
  mean_INS  <- mean(expr_matrix["INS", ])
  mean_NKX61 <- mean(expr_matrix["NKX6-1", ])

  # Fetch per-cell expression and cluster ID
  cell_expr <- FetchData(seurat_obj, vars = c("GCG", "ARX", "INS", "NKX6-1"))
  cell_expr$cluster <- seurat_obj[[cluster_col]][,1]

  # Calculate average expression per cluster
  cluster_avg <- aggregate(. ~ cluster, data = cell_expr, FUN = mean)

  # Identify clusters based on marker expression
  alpha_clusters <- cluster_avg$cluster[
    cluster_avg$GCG > mean_GCG &
      cluster_avg$ARX > mean_ARX &
      cluster_avg$INS < mean_INS ]

  beta_clusters <- cluster_avg$cluster[
    cluster_avg$INS > mean_INS &
      cluster_avg$`NKX6-1` > mean_NKX61 &
      cluster_avg$GCG < mean_GCG &
      !(cluster_avg$cluster %in% alpha_clusters) ]

  bihormonal_clusters <- cluster_avg$cluster[
    cluster_avg$GCG > mean_GCG &
      cluster_avg$INS > mean_INS &
      !(cluster_avg$cluster %in% c(alpha_clusters, beta_clusters)) ]

  # Assign cell types
  cluster_vector <- seurat_obj[[cluster_col]][,1]
  seurat_obj$CellType <- ifelse(
    cluster_vector %in% alpha_clusters, paste(prefix, "Alpha"),
    ifelse(
      cluster_vector %in% beta_clusters, paste(prefix, "Beta"),
      ifelse(
        cluster_vector %in% bihormonal_clusters, paste(prefix, "Bihormonal"),
        paste("Other", prefix)
      )
    )
  )

  return(seurat_obj)
}
