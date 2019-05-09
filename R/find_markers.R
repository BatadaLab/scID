#' @export
find_markers <- function(reference_gem, reference_clusters, logFC, only.pos=FALSE) {
  library(Seurat)
  
  so_ref <- CreateSeuratObject(reference_gem)
  so_ref <- suppressMessages(NormalizeData(so_ref))
  so_ref <- suppressMessages(ScaleData(so_ref))
  # For Seurat v2
  # so_ref@ident <- as.factor(reference_clusters)
  # For Seurat v3
  Idents(so_ref) <- as.factor(reference_clusters)
  markers <- suppressMessages(FindAllMarkers(so_ref, test.use = "MAST", only.pos = only.pos, logfc.threshold = logFC))
  
  return(markers)
}