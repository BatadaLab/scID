#' Main function to get Gene expression matrix and signature genes and return matches and scores
#' @param target_gem Data frame of gene expression (rows) per cell (columns) in target data
#' @param reference_gem Data frame of gene expression (rows) per cell (columns) in reference data
#' @param reference_clusters Named list of cluster IDs of reference cells
#' @param logFC LogFC threshold for extracting markers from reference clusters
#' @param use_reference_for_weights Logical to use either reference or target data for sorting the signature genes
#' @param likelihood_threshold Minimum required likelihood of gene signature score for a cell to be assigned to the respective reference cluster
#' @return list of target cluster IDs
#' @export
scid_match_cells <- function(target_gem = NULL, reference_gem = NULL, reference_clusters = NULL, 
                             logFC = 0.5, use_reference_for_weights = FALSE, likelihood_threshold = 0.99){
  
  #----------------------------------------------------------------------------------------------------
  # Check all reference cells have a cluster ID
  common_cells <- intersect(names(reference_clusters), colnames(reference_gem))
  if (length(common_cells) == 0) {
    print("None  of the reference cell has a cluster ID. Please check the reference_clusters list provided.")
    return()
  } else {
    reference_gem <- reference_gem[, common_cells]
    reference_clusters <- reference_clusters[common_cells]
  }
  # Find signature genes from reference data
  so_ref <- Seurat::CreateSeuratObject(raw.data = reference_gem)
  so_ref <- Seurat::NormalizeData(so_ref)
  so_ref <- Seurat::ScaleData(so_ref)
  so_ref@ident <- as.factor(reference_clusters)
  
  markers <- Seurat::FindAllMarkers(so_ref, test.use = "MAST", only.pos = TRUE, logfc.threshold = logFC)
  
  # Filter out signature genes that are not present in the target data
  markers <- markers[which(markers$gene %in% rownames(target_gem)), ]
  
  celltypes <- unique(markers$cluster)
  
  #----------------------------------------------------------------------------------------------------
  # Min-max normalization of target gem
  target_gem_norm <- t(apply(target_gem[markers$gene, ], 1, function(x) normalize_gem(x)))
  na.values <- which(apply(target_gem_norm, 1, function(x){any(is.na(x))}))
  if (length(na.values > 0)) {
    target_gem_norm <- target_gem_norm[-na.values, ]
  }
  
  #----------------------------------------------------------------------------------------------------
  # Weight signature genes
  weights <- list()
  if (use_reference_for_weights) {
    #Min-max normalization of reference gem
    ref_gem_norm <-  t(apply(reference_gem[markers$gene, ], 1, function(x) normalize_gem(x)))
    na.values <- which(apply(ref_gem_norm, 1, function(x){any(is.na(x))}))
    if (length(na.values > 0)) {
      ref_gem_norm <- ref_gem_norm[-na.values, ]
    }
    for (celltype in celltypes) {
      signature_genes <- markers$gene[which(markers$cluster == celltype)]
      IN <- names(which(reference_clusters == celltype))
      OUT <- setdiff(colnames(reference_gem), IN)
      gene.weights <- scID_weight(ref_gem_norm, IN, OUT)
      weights[[celltype]] <- gene.weights
    }
  } else {
    for (celltype in celltypes) {
      signature_genes <- markers$gene[which(markers$cluster == celltype)]
      putative_groups <- choose_unsupervised(target_gem, signature_genes)
      gene.weights <- scID_weight(target_gem_norm[signature_genes, ], putative_groups$in_pop, putative_groups$out_pop)
      weights[[celltype]] <- gene.weights
    }
  }
  
  #----------------------------------------------------------------------------------------------------
  # Find scores and putative matches
  scores <- data.frame(matrix(NA, length(celltypes), ncol(target_gem)), row.names = celltypes)
  colnames(scores) <- colnames(target_gem)
  for (celltype in celltypes) {
    signature <- names(weights[[celltype]])
    weighted_gem <- weights[[celltype]] * target_gem_norm[signature, ]
    score <- colSums(weighted_gem)/sum(weights[[celltype]])
    matches <- final_populations(score, likelihood_threshold)
    scores[celltype, matches] <- scale(score[matches])
  }
  
  # Resolve multiclass assignments
  labels <- c()
  for (cell in colnames(scores)) {
    if (all(is.na(scores[, cell]))) {
      matching_type <- "unassigned"
    } else {
      matching_type <- rownames(scores)[which(scores[, cell] == max(scores[, cell], na.rm = T))]
    }
    labels <- c(labels, matching_type)
  }
  names(labels) <- colnames(scores)
  table(labels)

  #----------------------------------------------------------------------------------------------------
  # return result
  return(labels)
    
  }
}
