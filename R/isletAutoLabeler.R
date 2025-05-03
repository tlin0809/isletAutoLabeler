#' Auto-label islet cell types
#'
#' Automatically labels Alpha, Beta, Bihormonal, Delta, PP, Duct, Acinar, Endothelial, Mesenchymal, Immune and Other cell clusters in a Seurat object using average marker expression.
#'
#' @param seurat_obj A Seurat object to be labeled.
#' @param assay The assay to use (default = "integrated")
#' @param cluster_col Column in metadata to use for clustering (default = "seurat_clusters")
#' @param source Cell source, used as a prefix for plot title (optional)
#' @param plot Logical (default = FALSE). If TRUE, generates and returns a labeled UMAP plot using DimPlot. To set plot = TRUE, there must exists a dimentionality reduction named umap. If FALSE, only returns Seurat object.
#' @return Updated Seurat object with a new column named "CellType", or Seurat object and labeled uMap if plot = TRUE.
#' @export

isletAutoLabeler <- function(seurat_obj, assay = "integrated",
                             cluster_col = "seurat_clusters", source = "Human",
                             plot = FALSE) {

  # Set prefix for plot title
  prefix <- source

  # Extract expression matrices
  expr_matrix <- GetAssayData(seurat_obj, assay = assay, slot = "data")

  # Calculate global average expression for each marker
  mean_GCG  <- mean(expr_matrix["GCG", ])
  mean_ARX  <- mean(expr_matrix["ARX", ])
  mean_INS  <- mean(expr_matrix["INS", ])
  mean_NKX61 <- mean(expr_matrix["NKX6-1", ])
  mean_SST <- mean(expr_matrix["SST", ])
  mean_PPY <- mean(expr_matrix["PPY", ])
  mean_KRT19 <- mean(expr_matrix["KRT19", ])
  mean_CPA1 <- mean(expr_matrix["CPA1", ])
  mean_PECAM1 <- mean(expr_matrix["PECAM1", ])
  mean_PDGFRB  <- mean(expr_matrix["PDGFRB", ])
  mean_FCER1G  <- mean(expr_matrix["FCER1G", ])

  # Fetch per-cell expression and cluster ID
  cell_expr <- FetchData(seurat_obj, vars = c("GCG", "ARX", "INS", "NKX6-1", "SST", "PPY", "KRT19", "CPA1", "PECAM1", "PDGFRB", "FCER1G"))
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

  delta_clusters <- cluster_avg$cluster[
    cluster_avg$SST > mean_SST &
      !(cluster_avg$cluster %in% c(alpha_clusters, beta_clusters, bihormonal_clusters)) ]

  pp_clusters <- cluster_avg$cluster[
    cluster_avg$PPY > mean_PPY &
      !(cluster_avg$cluster %in% c(alpha_clusters, beta_clusters, bihormonal_clusters, delta_clusters)) ]

  duct_clusters <- cluster_avg$cluster[
    cluster_avg$KRT19 > mean_KRT19 &
      !(cluster_avg$cluster %in% c(alpha_clusters, beta_clusters, bihormonal_clusters, delta_clusters, pp_clusters)) ]

  acinar_clusters <- cluster_avg$cluster[
    cluster_avg$CPA1 > mean_CPA1 &
      !(cluster_avg$cluster %in% c(alpha_clusters, beta_clusters, bihormonal_clusters, delta_clusters, pp_clusters, duct_clusters)) ]

  endothelial_clusters <- cluster_avg$cluster[
    cluster_avg$PECAM1 > mean_PECAM1 &
      !(cluster_avg$cluster %in% c(alpha_clusters, beta_clusters, bihormonal_clusters, delta_clusters, pp_clusters, duct_clusters, acinar_clusters)) ]

  mesenchymal_clusters <- cluster_avg$cluster[
    cluster_avg$PDGFRB > mean_PDGFRB &
      !(cluster_avg$cluster %in% c(alpha_clusters, beta_clusters, bihormonal_clusters, delta_clusters, pp_clusters, duct_clusters, acinar_clusters, endothelial_clusters)) ]

  immune_clusters <- cluster_avg$cluster[
    cluster_avg$FCER1G > mean_FCER1G &
      !(cluster_avg$cluster %in% c(alpha_clusters, beta_clusters, bihormonal_clusters, delta_clusters, pp_clusters, duct_clusters, acinar_clusters, endothelial_clusters, mesenchymal_clusters)) ]

  # Assign cell types
  cluster_vector <- seurat_obj[[cluster_col]][,1]
  seurat_obj$CellType <- ifelse(
    cluster_vector %in% alpha_clusters, "Alpha",
    ifelse(
      cluster_vector %in% beta_clusters, "Beta",
      ifelse(
        cluster_vector %in% bihormonal_clusters, "Bihormonal",
        ifelse(
          cluster_vector %in% delta_clusters, "Delta",
          ifelse(
            cluster_vector %in% pp_clusters, "PP",
            ifelse(
              cluster_vector %in% duct_clusters, "Duct",
              ifelse(
                cluster_vector %in% acinar_clusters, "Acinar",
                ifelse(
                  cluster_vector %in% endothelial_clusters, "Endothelial",
                  ifelse(
                    cluster_vector %in% mesenchymal_clusters, "Mesenchymal",
                    ifelse(
                      cluster_vector %in% immune_clusters, "Immune", "Other"
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )

  if (plot) {
    colors <- c(
      "Alpha" = "salmon",
      "Beta" = "deepskyblue",
      "Bihormonal" = "olivedrab3",
      "Delta" = "gold",
      "PP" = "plum",
      "Duct" = "slateblue2",
      "Acinar" = "bisque2",
      "Endothelial" = "slategray3",
      "Mesenchymal" = "yellow",
      "Immune" = "seagreen3",
      "Other" = "gray"
    )

  plot <- DimPlot(seurat_obj, group.by = "CellType", label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 5,
                  cols = colors) + theme(legend.position = 'none') +
    theme(legend.position = 'none',
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.9, "lines")) +
    ggtitle(paste(prefix, "-islet", sep = ""))

  return(list(object = seurat_obj, plot = plot))
  } else{
    return(seurat_obj)
  }

}
