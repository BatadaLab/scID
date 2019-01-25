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
  # Data pre-processing
  # Check all reference cells have a cluster ID
  common_cells <- intersect(names(reference_clusters), colnames(reference_gem))
  if (length(common_cells) == 0) {
    print("None  of the reference cells has a cluster ID. Please check the reference_clusters list provided.")
    return()
  } else {
    reference_gem <- reference_gem[, common_cells]
    reference_clusters <- reference_clusters[common_cells]
  }
  
  #----------------------------------------------------------------------------------------------------
  print("Stage 1: extract signatures genes from reference clusters")
  
  # Stage 1: Find signature genes from reference data
  so_ref <- Seurat::CreateSeuratObject(raw.data = reference_gem)
  so_ref <- Seurat::NormalizeData(so_ref)
  so_ref <- Seurat::ScaleData(so_ref)
  so_ref@ident <- as.factor(reference_clusters)
  
  markers <- suppressMessages(Seurat::FindAllMarkers(so_ref, test.use = "MAST", only.pos = TRUE, logfc.threshold = logFC))
  
  # Filter out signature genes that are not present in the target data
  markers <- markers[which(markers$gene %in% rownames(target_gem)), ]
  
  celltypes <- unique(markers$cluster)
  
  # Min-max normalization of target gem
  target_gem_norm <- t(apply(target_gem[markers$gene, ], 1, function(x) normalize_gem(x)))
  na.values <- which(apply(target_gem_norm, 1, function(x){any(is.na(x))}))
  if (length(na.values > 0)) {
    target_gem_norm <- target_gem_norm[-na.values, ]
  }
  
  #----------------------------------------------------------------------------------------------------
  # Stage 2: Weight signature genes
  print("Stage 2: Estimate weights of signature genes")
  
  weights <- list()
  
  if (use_reference_for_weights) {
    #Min-max normalization of reference gem
    ref_gem_norm <-  t(apply(reference_gem[markers$gene, ], 1, function(x) normalize_gem(x)))
    na.values <- which(apply(ref_gem_norm, 1, function(x){any(is.na(x))}))
    if (length(na.values > 0)) {
      ref_gem_norm <- ref_gem_norm[-na.values, ]
    }
    
    for (i in 1:length(celltypes)) {
      #svMisc::progress(i, max.value = length(celltypes), progress.bar = T)
      svMisc::progress(i*100/length(celltypes), max.value = 100, char = "%")
      Sys.sleep(0.01)
      signature_genes <- markers$gene[which(markers$cluster == celltypes[i])]
      IN <- names(which(reference_clusters == celltypes[i]))
      OUT <- setdiff(colnames(reference_gem), IN)
      gene.weights <- scID_weight(ref_gem_norm, IN, OUT)
      weights[[celltypes[i]]] <- gene.weights
      if (i==length(celltypes)) cat("Done!")
    }
    
  } else {
    for (i in 1:length(celltypes)) {
      svMisc::progress(i*100/length(celltypes), max.value = 100, char = "%")
      Sys.sleep(0.01)
      signature_genes <- markers$gene[which(markers$cluster == celltypes[i])]
      putative_groups <- choose_unsupervised(target_gem, signature_genes)
      gene.weights <- scID_weight(target_gem_norm[signature_genes, ], putative_groups$in_pop, putative_groups$out_pop)
      weights[[celltypes[i]]] <- gene.weights
      if (i==length(celltypes)) cat("Done!")
    }
  }
  
  #----------------------------------------------------------------------------------------------------
  # Stage 3
  # Find scores and putative matches
  print("Stage 3.1-2: Calculate scores and find matching cells")
  
  scores <- data.frame(matrix(NA, length(celltypes), ncol(target_gem)), row.names = celltypes)
  colnames(scores) <- colnames(target_gem)
  
  for (i in 1:length(celltypes)) {
    celltype <- as.integer(celltypes[i])
    svMisc::progress(i*100/length(celltypes), max.value = 100)
    Sys.sleep(0.01)
    signature <- names(weights[[celltype]])
    weighted_gem <- weights[[celltype]] * target_gem_norm[signature, ]
    score <- colSums(weighted_gem)/sum(weights[[celltype]])
    matches <- final_populations(score, likelihood_threshold) # NEED TO SUPPRESS THE MESSAGE
    scores[celltype, matches] <- scale(score[matches])
    if (i==length(celltypes)) cat("Done!")
  }
  
  # Resolve multiclass assignments
  print ("Stage 3.3: Resolve multiclass assignments")
  labels <- c()
  for (i in 1:ncol(scores)) {
    svMisc::progress(i*100/ncol(scores), max.value = 100)
    Sys.sleep(0.01)
    cell <- colnames(scores)[i]
    if (all(is.na(scores[, cell]))) {
      matching_type <- "unassigned"
    } else {
      matching_type <- rownames(scores)[which(scores[, cell] == max(scores[, cell], na.rm = T))]
    }
    labels <- c(labels, matching_type)
    if (i==ncol(scores)) cat("Done!")
  }
  names(labels) <- colnames(scores)
  table(labels)

  # return result
  return(labels)
    
}
