#' @export
find_markers <- function(reference_gem, reference_clusters, logFC) {
  so_ref <- Seurat::CreateSeuratObject(raw.data = reference_gem)
  so_ref <- suppressMessages(Seurat::NormalizeData(so_ref))
  so_ref <- suppressMessages(Seurat::ScaleData(so_ref))
  so_ref@ident <- as.factor(reference_clusters)
  
  markers <- suppressMessages(Seurat::FindAllMarkers(so_ref, test.use = "MAST", only.pos = TRUE, logfc.threshold = logFC))
  
  return(markers)
}