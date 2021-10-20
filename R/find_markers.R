#' Function to extract cluster specific genes from reference clusters.
#' This function uses MAST test as implemented in the Seurat package.
#' 
#' @param reference_gem Data frame of gene expression (rows) per cell (columns) in reference data
#' @param reference_clusters Named list of cluster IDs of reference cells
#' @param logFC Log2FC threshold for extracting markers from reference clusters
#' @param only.pos Logical to include negative markers in the cluster specific gene sets
#' @param normalize_reference Logical to select if reference data need to be normalized (when raw counts have been provided)
#' 
#' @return Data frame of cluster specific genes extracted
#' 
#' @export
find_markers <- function(reference_gem, reference_clusters, logFC, only.pos, normalize_reference) {
  library(Seurat)
  
  so_ref <- CreateSeuratObject(reference_gem)
  if (normalize_reference) {
    so_ref <- suppressMessages(NormalizeData(so_ref))
  }
  so_ref <- suppressMessages(ScaleData(so_ref))
  Idents(so_ref) <- as.factor(reference_clusters)
  markers <- suppressMessages(FindAllMarkers(so_ref, test.use = "MAST", only.pos = only.pos, logfc.threshold = logFC))
  
  markers
}
