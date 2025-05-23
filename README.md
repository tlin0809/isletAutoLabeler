# One-line Labeling islet cell-types
Automatically labels Alpha, Beta, Bihormonal, Delta, PP, Duct, Acinar, Endothelial, Mesenchymal, Immune, and Other cell clusters in a Seurat object using average marker expression. Provides a plotting option for quick visualization of the labeled uMap of the updated object. 


# Install dependencies (if not already installed)
install.packages("remotes")

# Install the isletAutoLabeler package
remotes::install_github("tlin0809/isletAutoLabeler")

# Load isletAutoLabeler: 
library(isletAutoLabeler)

# Load help page: 
?isletAutoLabeler

# Usage: 
isletAutoLabeler(
          seurat_obj,
          assay = "integrated",
          cluster_col = "seurat_clusters",
          source = "Human",
          plot = FALSE
)

# Notes:
Use isletAutoLabeler() after dimensionality reduction (e.g., PCA), clustering, and running RunUMAP(). 
isletAutoLabeler() assumes the passed Seurat object is already in label-ready form.


