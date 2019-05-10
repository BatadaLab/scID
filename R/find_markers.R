#' @export
find_markers <- function(reference_gem, reference_clusters, logFC) {
  library(Seurat)
  so_ref <- CreateSeuratObject(counts = reference_gem)
  so_ref <- suppressMessages(NormalizeData(so_ref))
  so_ref <- suppressMessages(ScaleData(so_ref))
  Idents(so_ref) <- as.factor(reference_clusters)
  
  markers <- suppressMessages(FindAllMarkers(so_ref, test.use = "MAST", only.pos = TRUE, logfc.threshold = logFC))
  
  return(markers)
}
