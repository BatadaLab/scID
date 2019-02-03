#' @export
make_heatmap <- function(gem, labels, markers) {
  
  rownames(gem) <- toupper(rownames(gem))

  # Keep markers present in gem
  markers <- markers[which(markers$gene %in% rownames(gem)), ]
  
  celltypes <- unique(labels)
  
  nrows <- length(unique(markers$cluster))
  ncols <- length(celltypes)
  
  gem_avg <- matrix(NA, nrows, ncols)
  for (i in 1:ncols) {
    cells <- na.omit(names(labels)[which(labels == celltypes[i])])
    if (length(cells) > 1) {
      avg_exp <- rowMeans(gem[toupper(markers$gene), cells])
    } else if (length(cells) == 1) {
      avg_exp <- gem[toupper(markers$gene), cells]
      names(avg_exp) <- toupper(markers$gene)
    } else {
      next
    }
    for (j in 1:nrows) {
      gem_avg[j,i] <- mean(na.omit(avg_exp[toupper(markers$gene[which(markers$cluster == unique(markers$cluster)[j])])]))
    }
  }
  
  rownames(gem_avg) <- paste("gs", unique(markers$cluster), sep = "_")
  colnames(gem_avg) <- paste("Cl", celltypes, sep = "_")
  
  pheatmap::pheatmap(gem_avg, border="white", #color = colorspace::diverge_hcl(50, h=c(180, 70), c=70, l=c(70, 90)), 
                     cluster_rows = F, cluster_cols = F, border_color = F, scale = "row")
}