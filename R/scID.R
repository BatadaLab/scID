#' Main function to get Gene expression matrix and signature genes and return matches and scores
#' @param target_gem Data frame of gene expression (rows) per cell (columns) in target data
#' @param reference_gem Data frame of gene expression (rows) per cell (columns) in reference data
#' @param reference_clusters Named list of cluster IDs of reference cells
#' @param logFC LogFC threshold for extracting markers from reference clusters
#' @param use_reference_for_weights Logical to use either reference or target data for sorting the signature genes
#' @param likelihood_threshold Minimum required likelihood of gene signature score for a cell to be assigned to the respective reference cluster
#' @return list of target cluster IDs
#' @export
scid_match_cells <- function(target_gem = NULL, reference_gem = NULL, reference_clusters = NULL, markers = NULL,
                             logFC = 0.5, estimate_weights = TRUE, weights = NULL, only_pos=FALSE){
  
  #----------------------------------------------------------------------------------------------------
  # Data pre-processing
  if (is.null(reference_gem) && is.null(reference_clusters) && is.null(markers)) {
    message("Please provide either clustered reference data or list of markers for each reference cluster")
    return()
  } else if (is.null(markers)) {
    # Check all reference cells have a cluster ID
    common_cells <- intersect(names(reference_clusters), colnames(reference_gem))
    if (length(common_cells) == 0) {
      print("None  of the reference cells has a cluster ID. Please check the reference_clusters list provided.")
      return()
    } else {
      reference_gem <- reference_gem[, common_cells]
      rownames(reference_gem) <- toupper(rownames(reference_gem))
      # Remove genes that are zero across all cells
      reference_gem <- reference_gem[which(rowSums(reference_gem) != 0), ]
      reference_clusters <- reference_clusters[common_cells]
    }
  } else {
    # Check markers have gene and cluster columns
    if (length(intersect(c("gene", "cluster"), colnames(markers))) !=2 ) {
      message("Please provide a data frame of markers with gene and cluster in columns")
      return() 
    }
    if (!"gene" %in% colnames(markers)) {
      message("Please provide a data frame of markers with gene and cluster in columns")
      return()  
    } else {
      markers$gene <- toupper(markers$gene)
    }
  }

  # Target
  rownames(target_gem) <- toupper(rownames(target_gem))
  # Remove genes that are zero across all cells
  target_gem <- target_gem[which(rowSums(target_gem) != 0), ]

  # ----------------------------------------------------------------------------------------------------
  # Stage 1: Find signature genes from reference data
  if (is.null(markers)) {
    # ----------------------------------------------------------------------------------------------------
    # Stage 1: Find signature genes from reference data
    message("Stage 1: extract signatures genes from reference clusters")
    markers <- find_markers(reference_gem, reference_clusters, logFC, only.pos=only_pos)
    # Filter out signature genes that are not present in the target data
    markers <- markers[which(markers$gene %in% rownames(target_gem)), ]
    celltypes <- unique(markers$cluster)
  } else {
    markers <- markers[which(markers$gene %in% rownames(target_gem)), ]
    celltypes <- unique(markers$cluster)
  }
  
  # # ----------------------------------------------------------------------------------------------------
  # # Assess markers with self-mapping
  # print("Assessing quality of extracted markers")
  # # ref_cells_sbst <- c()
  # # for (celltype in celltypes) {
  # #   ref_cells_sbst <- c(ref_cells_sbst, sample(Seurat::WhichCells(so_ref, celltype), round(0.2 * length(Seurat::WhichCells(so_ref, celltype)))))
  # # }
  # # ref_subset <- reference_gem[, ref_cells_sbst]
  # 
  # weights <- list()
  # # Min-max normalization of reference gem
  # ref_gem_norm <-  t(apply(reference_gem[markers$gene, ], 1, function(x) normalize_gem(x)))
  # na.values <- which(apply(ref_gem_norm, 1, function(x){any(is.na(x))}))
  # if (length(na.values > 0)) {
  #   ref_gem_norm <- ref_gem_norm[-na.values, ]
  # }
  # 
  # for (i in 1:length(celltypes)) {
  #   svMisc::progress(i*100/length(celltypes), max.value = 100, char = "-", progress.bar = T)
  #   Sys.sleep(0.01)
  #   signature_genes <- markers$gene[which(markers$cluster == celltypes[i])]
  #   IN <- intersect(names(which(reference_clusters == celltypes[i])), colnames(ref_gem_norm))
  #   OUT <- setdiff(colnames(ref_gem_norm), IN)
  #   gene.weights <- scID_weight(ref_gem_norm, IN, OUT)
  #   weights[[celltypes[i]]] <- gene.weights
  #   if (i==length(celltypes)) cat("Done!")
  # }
  # 
  # scores <- data.frame(matrix(NA, length(celltypes), ncol(ref_gem_norm)), row.names = celltypes)
  # colnames(scores) <- colnames(ref_gem_norm)
  # for (i in 1:length(celltypes)) {
  #   celltype <- celltypes[i]
  #   svMisc::progress(i*100/length(celltypes), max.value = 100, char = "-", progress.bar = T)
  #   Sys.sleep(0.01)
  #   signature <- names(weights[[celltype]])
  #   weighted_gem <- weights[[celltype]] * ref_gem_norm[signature, ]
  #   score <- colSums(weighted_gem)/sum(weights[[celltype]])
  #   matches <- final_populations(score, likelihood_threshold) # NEED TO SUPPRESS THE MESSAGE
  #   scores[as.character(celltype), matches] <- scale(score[matches])
  #   if (i==length(celltypes)) cat("Done!")
  # }
  # 
  # labels <- c()
  # for (i in 1:ncol(scores)) {
  #   svMisc::progress(i*100/ncol(scores), max.value = 100, char = "-", progress.bar = T)
  #   Sys.sleep(0.01)
  #   cell <- colnames(scores)[i]
  #   if (all(is.na(scores[, cell]))) {
  #     matching_type <- "unassigned"
  #   } else {
  #     matching_type <- rownames(scores)[which(scores[, cell] == max(scores[, cell], na.rm = T))]
  #   }
  #   labels <- c(labels, matching_type)
  #   if (i==ncol(scores)) cat("Done!")
  # }
  # names(labels) <- colnames(scores)
  # 
  # ARI <- mclust::adjustedRandIndex(labels, so_ref@ident[names(labels)])
  # 
  # if (ARI < 0.5) {
  #   print("scID could not extract discriminative markers from reference clusters provided. Please try again with different clusters or logFC threshold.")
  #   return()
  # } else {
  #   print(paste("ARI of self-mapping is ", ARI, sep = ""))
  # }

  # ----------------------------------------------------------------------------------------------------
  # Min-max normalization of target gem
  target_gem_norm <- t(apply(target_gem[markers$gene, ], 1, function(x) normalize_gem(x)))
  na.values <- which(apply(target_gem_norm, 1, function(x){any(is.na(x))}))
  if (length(na.values > 0)) {
    target_gem_norm <- target_gem_norm[-na.values, ]
  }
  
  # ----------------------------------------------------------------------------------------------------
  # Stage 2: Weight signature genes
  if (estimate_weights) {
    print("Stage 2: Estimate weights of signature genes")
    weights <- list()
    training_labels <- data.frame(matrix(NA, nrow = ncol(target_gem), ncol = length(celltypes)), row.names = colnames(target_gem))
    colnames(training_labels) <- celltypes
    
    for (i in 1:length(celltypes)) {
      celltype_markers <- markers[which(markers$cluster == celltypes[i]), ]
      positive_markers <- celltype_markers$gene[which(celltype_markers$avg_logFC > 0)]
      negative_markers <- celltype_markers$gene[which(celltype_markers$avg_logFC < 0)]
      
      putative_groups <- choose_unsupervised(target_gem, positive_markers, negative_markers)
      training_labels[putative_groups$in_pop, i] <- "IN"
      training_labels[putative_groups$out_pop, i] <- "OUT"
      signature_genes <- c(positive_markers, negative_markers)
      gene.weights <- scID_weight(target_gem_norm[signature_genes, ], putative_groups$in_pop, putative_groups$out_pop)
      weights[[as.character(celltypes[i])]] <- gene.weights
      svMisc::progress(i*100/length(celltypes), max.value = 100, char = "-", progress.bar = T)
      Sys.sleep(0.01)
      if (i==length(celltypes)) cat("Done!")
    }
    names(weights) <- celltypes
  } else {
    if (is.null(weights)) {
      message("Please provide weights of signature genes or choose to estimate weights from target data")
      return()
    }
    # for (i in 1:length(celltypes)) {
    #   signature_genes <- markers$gene[which(markers$cluster == celltypes[i])]
    #   gene.weights <- markers$weight[which(markers$cluster == celltypes[i])]
    #   names(gene.weights) <- signature_genes
    #   weights[[celltypes[i]]] <- gene.weights
    # }
  }
  
  #----------------------------------------------------------------------------------------------------
  # Stage 3
  # Find scores and putative matches
  print("Stage 3.1-2: Calculate scores and find matching cells")
  
  scores <- data.frame(matrix(NA, length(celltypes), ncol(target_gem)), row.names = celltypes)
  colnames(scores) <- colnames(target_gem)
  
  for (i in 1:length(celltypes)) {
    celltype <- as.character(celltypes[i])
    svMisc::progress(i*100/length(celltypes), max.value = 100, char = "-", progress.bar = T)
    Sys.sleep(0.01)
    signature <- intersect(names(weights[[celltype]]), rownames(target_gem_norm))
    weighted_gem <- weights[[celltype]][signature] * target_gem_norm[signature, ]
    score <- colSums(weighted_gem)/sqrt(sum(weights[[celltype]]^2))
    scores[as.character(celltype), ] <- score
    #matches <- final_populations(score) 
    #scores[as.character(celltype), matches] <- scale(score[matches])
    if (i==length(celltypes)) cat("Done!")
  }
  return(list(scores=scores, markers=markers, weights=weights))
  
  # Resolve multiclass assignments
  print ("Stage 3.3: Resolve multiclass assignments")
  labels <- c()
  for (i in 1:ncol(scores)) {
    svMisc::progress(i*100/ncol(scores), max.value = 100, char = "-", progress.bar = T)
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

  # return result
  return(list(scores=scores, markers=markers, estimated_weights=weights, training_labels=training_labels, labels=labels))

}
