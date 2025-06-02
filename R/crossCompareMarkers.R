#' Cross compare marker gene lists across two comparisons
#'
#' @title crossCompareMarkers
#' @description Filters and cross-compares two sets of potential markers derived from findMechanismMarkers() and
#' identifies potential markers with shared and significant directionality in expression across comparisons.
#'
#' @param result1_scg Output strongly-correlated gene list from findMechanismMarkers(), comparison 1 (specify raw or filtered).
#' @param result2_scg Output strongly-correlated gene list from findMechanismMarkers(), comparison 2 (specify raw or filtered).
#' @param deg_comp1 Data frame. Full DEG table from findMechanismMarkers(), comparison 1. Used to compare expression direction of genes from result2.
#' @param deg_comp2 Data frame. Full DEG table from findMechanismMarkers(), comparison 2. Used to compare expression direction of genes from result1.
#' @param direction_shared Character. Either "down" or "up" for filtering by shared expression pattern (default = down).
#' @param logfc_threshold Numeric. Log2 fold change threshold for defining directionality (default = 0).
#' @param sig_threshold Numeric. p-value threshold for keeping significant DEGs (default = 0.05).
#'
#' @return A list containing:
#'   - markers_1: character vector of candidate marker genes from comparison 1
#'   - markers_2: character vector of candidate marker genes from comparison 2
#'   - shared_markers: intersected genes from the two potential marker lists
#'   - scg1_table: Data frame tracking filtering of comparison 1's SCGs through comp1 and comp2
#'   - scg2_table: Data frame tracking filtering of comparison 2's SCGs through comp2 and comp1
#' @export

crossCompareMarkers <- function(
    result1_scg,
    result2_scg,
    deg_comp1,
    deg_comp2,
    direction_shared = "down",
    logfc_threshold = 0,
    sig_threshold = 0.05
) {
  # Extract strongly correlated gene lists
  scg1 <- unique(unlist(result1_scg))
  scg2 <- unique(unlist(result2_scg))

  # Confirm scg1 expression in comparison 1
  scg1_comp1 <- deg_comp1[rownames(deg_comp1) %in% scg1, ]
  if (direction_shared == "down") {
    scg1_comp1 <- subset(scg1_comp1, avg_log2FC < logfc_threshold & p_val < sig_threshold)
  } else if (direction_shared == "up") {
    scg1_comp1 <- subset(scg1_comp1, avg_log2FC > logfc_threshold & p_val < sig_threshold)
  }

  # Then check scg1_comp1 expression in comparison 2
  scg1_comp2 <- deg_comp2[rownames(deg_comp2) %in% rownames(scg1_comp1), ]
  if (direction_shared == "down") {
    scg1_comp2 <- subset(scg1_comp2, avg_log2FC < logfc_threshold & p_val < sig_threshold)
  } else if (direction_shared == "up") {
    scg1_comp2 <- subset(scg1_comp2, avg_log2FC > logfc_threshold & p_val < sig_threshold)
  }
  potential_markers_1 <- rownames(scg1_comp2)

  # Confirm scg2 expression in comparison 2
  scg2_comp2 <- deg_comp2[rownames(deg_comp2) %in% scg2, ]
  if (direction_shared == "down") {
    scg2_comp2 <- subset(scg2_comp2, avg_log2FC < logfc_threshold & p_val < sig_threshold)
  } else if (direction_shared == "up") {
    scg2_comp2 <- subset(scg2_comp2, avg_log2FC > logfc_threshold & p_val < sig_threshold)
  }

  # Then check scg2_comp2 expression in comparison 1
  scg2_comp1 <- deg_comp1[rownames(deg_comp1) %in% rownames(scg2_comp2), ]
  if (direction_shared == "down") {
    scg2_comp1 <- subset(scg2_comp1, avg_log2FC < logfc_threshold & p_val < sig_threshold)
  } else if (direction_shared == "up") {
    scg2_comp1 <- subset(scg2_comp1, avg_log2FC > logfc_threshold & p_val < sig_threshold)
  }
  potential_markers_2 <- rownames(scg2_comp1)

  # Find shared markers
  shared_potential_markers <- intersect(potential_markers_1, potential_markers_2)

  # Summary tables for tracing filtering steps
  scg1_table <- data.frame(
    orig_scg1 = scg1,
    retained_in_comp1 = scg1 %in% rownames(scg1_comp1),
    retained_in_comp2 = scg1 %in% rownames(scg1_comp2)
  )

  scg2_table <- data.frame(
    orig_scg2 = scg2,
    retained_in_comp2 = scg2 %in% rownames(scg2_comp2),
    retained_in_comp1 = scg2 %in% rownames(scg2_comp1)
  )

  return(list(
    markers_1 = potential_markers_1,
    markers_2 = potential_markers_2,
    shared_markers = shared_potential_markers,
    scg1_table = scg1_table,
    scg2_table = scg2_table
  ))
}
