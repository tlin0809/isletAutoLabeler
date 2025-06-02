#' @title crossCompareMarkers
#' @description This function takes the output of two findMechanismMarkers() runs and identifies shared genes that are downregulated and strongly correlated in both comparisons.
#' @param result1 Output list from findMechanismMarkers() for comparison 1.
#' @param result2 Output list from findMechanismMarkers() for comparison 2.
#' @param marker_genes Character vector. Marker genes used in both comparisons.
#' @param corr_threshold Numeric. Absolute R threshold for selecting strongly correlated genes. Default is 0.5.
#' @return A list containing:
#'   - potential_markers_1: character vector of candidate marker genes from comparison 1
#'   - potential_markers_2: character vector of candidate marker genes from comparison 2
#'   - shared_potential_markers: intersected genes from the two potential marker lists
#' @export
crossCompareMarkers <- function(result1, result2, marker_genes, corr_threshold = 0.5) {
  # Extract gene names strongly correlated with each marker
  strong_1 <- unique(unlist(result1$strongly_correlated_genes[marker_genes]))
  strong_2 <- unique(unlist(result2$strongly_correlated_genes[marker_genes]))

  # Extract gene names that are in the downregulated DEGs
  de_genes_1 <- rownames(subset(result1$top_markers, avg_log2FC < 0))
  de_genes_2 <- rownames(subset(result2$top_markers, avg_log2FC < 0))

  # Get genes that are both strongly correlated AND downregulated
  potential_markers_1 <- intersect(strong_1, de_genes_1)
  potential_markers_2 <- intersect(strong_2, de_genes_2)

  # Intersect both sets to get shared potential markers
  shared_potential_markers <- intersect(potential_markers_1, potential_markers_2)

  return(list(
    potential_markers_1 = potential_markers_1,
    potential_markers_2 = potential_markers_2,
    shared_potential_markers = shared_potential_markers
  ))
}
