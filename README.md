# isletToolkit

`isletToolkit` is a modular R package designed for the single-cell RNA-seq analysis of pancreatic islet cells.  
It provides tools for:

- Automated cell type labeling using average marker gene expression
- Identification of potential mechanism-associated marker genes (e.g., maturation, proliferation) through differential expression and correlation analysis
- Cross-comparison of marker genes across developmental or experimental conditions to discover promising machanism-associated markers

## Installation

```r
# Install dependencies if not already installed
install.packages("remotes")

# Install the isletToolkit package from GitHub
remotes::install_github("tlin0809/isletToolkit")
```

## Loading Package and Help Page
```r
# Load isletAutoLabeler: 
library(isletAutoLabeler)
```

```r
# Load help page: 
?isletAutoLabeler
```

# Function Overview  

## isletAutoLabeler()
Automatically labels major islet cell types (Alpha, Beta, Bihormonal, Delta, PP, 
Duct, Acinar, Endothelial, Mesenchymal, Immune, and Other) in a Seurat object 
based on average marker gene expression. 

Optionally generates a UMAP plot of labeled clusters.

```r
isletAutoLabeler(
  seurat_obj,
  assay = "integrated",
  cluster_col = "seurat_clusters",
  source = "Human",
  plot = TRUE
)
```
Outputs include: 
- Updated Seurat object with added "CellType" column
- Labeled UMAP if plot = "TRUE"

## findMechanismMarkers()
Identifies candidate marker genes for a user-defined cellular mechanism (e.g., 
maturation, proliferation) by integrating differential expression (DE) and Pearson
correlation coefficient (PCC) analysis.
Genes are correlated with user-specified marker genes to identify candidates that
are both statistically and biologically relevant.

To run subsequent cross-comparision analysis, run findMechanismMarkers() twice 
with separate population groups.
```r
result <- findMechanismMarkers(
  seurat_obj,
  ident.1 = "SC-bihormonal",
  ident.2 = "SC-alpha",
  top_n_DEGs_volcano = 20,
  top_n_DEGs_pcc = 500,
  direction = "down",
  logfc.threshold = -0.5,
  pval.threshold = 0.05,
  marker_genes = c("ARX", "GCG", "TTR", "FOX2", "MAFB"),
  corr_threshold = 0.5
)
```
Outputs include: 
- Volcano plot of DEGs

- Filtered and unfiltered DEGs

- PCC matrices

- Lists of strongly correlated genes

- Correlation tables (raw and filtered)

## crossCompareMarkers()
Cross-compares the outputs from multiple findMechanismMarkers() calls to identify
marker genes consistently regulated across two developmental stages or treatment conditions.

This function increases biological confidence by ensuring that selected genes 
follow a consistent trend. For example, a maturation marker should show increasing
expression from SC-bihormonal to SC-alpha to Human-alpha.
```r
cross_results <- crossCompareMarkers(
  result1 = res1$strongly_correlated_genes_filtered,
  result2 = res2$strongly_correlated_genes_filtered,
  marker_genes = c("ARX", "GCG", "TTR", "FOX2", "MAFB"),
  direction_shared = "down",
  corr_threshold = 0.5
)
```
Outputs include:

- Marker gene list from comparison 1

- Marker gene list from comparison 2

- Intersected gene list between marker gene list 1 and marker gene list 2

# Example Workflow 
```r
# Step 1: Annotate cell types
seurat_obj <- isletAutoLabeler(seurat_obj)

# Step 2: Run mechanism marker discovery on two comparisons
res1 <- findMechanismMarkers(seurat_obj, ident.1 = "SC-bihormonal", ident.2 = "SC-alpha", ...)
res2 <- findMechanismMarkers(seurat_obj, ident.1 = "SC-alpha", ident.2 = "Human-alpha", ...)

# Step 3: Cross-compare markers from two developmental transitions
cross_results <- crossCompareMarkers(
  res1$strongly_correlated_genes_filtered,
  res2$strongly_correlated_genes_filtered,
  direction = "down"
)
```

# Complete Example

## Load the package
```r
library(isletToolkit)
```

## Step 1: Annotate cell types
```r
seurat_obj <- isletAutoLabeler(
  seurat_obj,
  assay = "integrated",
  cluster_col = "seurat_clusters",
  source = "Human",
  plot = TRUE
)
```

## Step 2: Run mechanism marker discovery on two (or more) comparisons
```r
res1 <- findMechanismMarkers(
  seurat_obj,
  ident.1 = "SC-bihormonal",
  ident.2 = "SC-alpha",
  top_n_DEGs_volcano = 20,
  top_n_DEGs_pcc = 500,
  direction = "down",
  logfc.threshold = -0.5,
  pval.threshold = 0.05,
  marker_genes = c("ARX", "GCG", "TTR", "FOX2", "MAFB"),
  corr_threshold = 0.5
)
```

```r
res2 <- findMechanismMarkers(
  seurat_obj,
  ident.1 = "SC-alpha",
  ident.2 = "Human-alpha",
  top_n_DEGs_volcano = 20,
  top_n_DEGs_pcc = 500,
  direction = "down",
  logfc.threshold = -0.5,
  marker_genes = c("ARX", "GCG", "TTR", "FOX2", "MAFB"),
  corr_threshold = 0.5
)
```

## Step 3: Cross-compare markers from two developmental transitions
```r
cross_results <- crossCompareMarkers(
  result1 = res1$strongly_correlated_genes_filtered,
  result2 = res2$strongly_correlated_genes_filtered,
  marker_genes = c("ARX", "GCG", "TTR", "FOX2", "MAFB"),
  direction_shared = "down",
  corr_threshold = 0.5
)
```

## Notes

- Use `isletAutoLabeler()` after dimensionality reduction (e.g., PCA), clustering, and running `RunUMAP()` on your Seurat object.
- The function assumes the object is label-ready and uses default human markers unless specified otherwise.
- The output of `findMechanismMarkers()` includes volcano plots, DEGs, correlation matrices, and filtered lists of strongly correlated candidate markers.
- `crossCompareMarkers()` helps reinforce biological relevance by filtering for direction and statistical significance across comparisons.

## License

MIT License

## Citation and Contributions

If you use this package in your work, please cite the GitHub repository.  
Contributions are welcome through pull requests or GitHub Issues.
