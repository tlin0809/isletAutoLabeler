Automatically labels Alpha, Beta, Bihormonal, Delta, PP, Duct, Acinar, Endothelial, Mesenchymal, Immune, and Other cell clusters in a Seurat object using average marker expression. Provides a plotting option for quick visualization of the labeled uMap of the updated object. 

Download isletAutoLabeler in any R session:  devtools::install_github("tlin0809/isletAutoLabeler")
Load isletAutoLabeler: library(isletAutoLabeler)
Load help page: ?isletAutoLabeler

Usage: isletAutoLabeler(
          seurat_obj,
          assay = "integrated",
          cluster_col = "seurat_clusters",
          source = "Human",
          plot = FALSE
        )


