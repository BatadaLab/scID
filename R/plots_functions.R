#' Function to plot heatmap of average cluster-specific geneset 
#' expression in clusters of cells
#' 
#' @param gem Data frame of gene expression (rows) per cell (columns) in target data
#' @param labels List of cluster IDs for cells of gem
#' @param markers Data frame of cluster specific genes with at least a "gene" and a "cluster" column
#' 
#' @export
make_heatmap <- function(gem, labels, markers) {
  
  # Keep only cells with available labels
  common_cells <- intersect(colnames(gem), names(labels))
  if (length(common_cells) == 0) {
    stop("Cell names between labels and gem do not match! Please make sure you have provided labels for the cells of the gem.")
  } else {
    gem <- gem[, common_cells]
    labels <- labels[common_cells]
  }
  
  rownames(gem) <- make.names(toupper(rownames(gem)), unique=TRUE)
  markers$gene <- toupper(markers$gene)

  # Keep positive markers
  markers <- markers[which(markers$avg_log2FC > 0), ]
  # Keep markers present in gem
  markers <- markers[which(markers$gene %in% rownames(gem)), ]
  
  celltypes <- unique(c(unique(as.character(markers$cluster)), unique(as.character(labels))))
  
  gem_avg <- matrix(NA, length(celltypes), length(celltypes))
  for (i in 1:length(celltypes)) {
    cells <- na.omit(names(labels)[which(labels == celltypes[i])])
    if (length(cells) >= 1) {
      avg_exp <- rowMeans(gem[markers$gene, cells, drop = FALSE])
    } else {
      next
    }
    for (j in 1:length(celltypes)) {
      gem_avg[j,i] <- mean(avg_exp[markers$gene[which(markers$cluster == celltypes[j])]], na.rm = TRUE)
    }
  }
  
  rownames(gem_avg) <- paste("gs", celltypes, sep = "_")
  colnames(gem_avg) <- paste("Cl", celltypes, sep = "_")
  # remove columns that are all NA
  na.cols = which(apply(gem_avg, 2, function(x) all(is.na(x))))
  if (length(na.cols) > 0) {
    gem_avg <- gem_avg[, -na.cols]
  }
  na.rows = which(apply(gem_avg, 1, function(x) all(is.na(x))))
  if (length(na.rows) > 0) {
    # remove rows that are all NA
    gem_avg <- gem_avg[-na.rows, ]
  }
  
  pheatmap::pheatmap(gem_avg, border="white", cluster_rows = F, cluster_cols = F, border_color = F, scale = "row")
}

#' Function to plot heatmap of average cluster-specific geneset 
#' expression in clusters of cells
#' 
#' @param gem Data frame of gene expression (rows) per cell (columns) in target data
#' @param labels List of cluster IDs for cells of gem
#' @param markers Data frame of cluster specific genes with at least "gene", "cluster" and "avg_logFC" columns
#' @param clusterID cluster ID of cluster of interest
#' @param weights list of weighst of cluster specific genes as returned by scid_multiclass
#' 
#' @export
plot_score_2D <- function(gem, labels, markers, clusterID, weights) {
  
  rownames(gem) <- make.names(toupper(rownames(gem)), unique=TRUE)
  markers$gene <- toupper(markers$gene)
  
  markers <- markers[which(markers$cluster == clusterID), ]
  positive_markers <- intersect(markers$gene[which(markers$avg_log2FC > 0)], rownames(gem))
  negative_markers <- intersect(markers$gene[which(markers$avg_log2FC < 0)], rownames(gem))
  
  if (length(positive_markers) == 0) {
    stop("No positive markers available for this cell type")
  } else if (length(negative_markers) == 0) {
    stop("No negative markers available for this cell type")
  } else {
    gem_norm <- t(apply(gem[c(positive_markers, negative_markers), ], 1, function(x) normalize_gene(x)))
    gem_norm <- gem_norm[complete.cases(gem_norm), ]
    
    weighted_gem <- weights[[clusterID]][c(positive_markers, negative_markers)] * gem_norm[, ,drop=FALSE]
    
    df <- data.frame(positive_score = colSums(weighted_gem[positive_markers, , drop=FALSE])/sqrt(sum(weights[[clusterID]][positive_markers]^2)), 
                     negative_score = colSums(weighted_gem[negative_markers, , drop=FALSE])/sqrt(sum(weights[[clusterID]][negative_markers]^2)))
    df$label <- rep("Other cell type", nrow(df))
    df[names(labels)[which(labels == clusterID)], "label"] <- clusterID
    df$label <- factor(df$label, levels = c(clusterID, "Other cell type"))
    
    library(ggplot2)
    ggplot(df, aes(x=positive_score, y=negative_score, color=label)) + geom_point() + 
      scale_color_manual(values=c("black", "grey")) + theme_classic()
  }
}


  
  
  
