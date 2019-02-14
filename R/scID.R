#' Main function to get Gene expression matrix and signature genes and return matches and scores
#' @param target_gem Data frame of gene expression (rows) per cell (columns) in target data
#' @param reference_gem (optional) Data frame of gene expression (rows) per cell (columns) in reference data
#' @param reference_clusters (optional) Named list of cluster IDs of reference cells
#' @param markers (optional) Data frame with of genes and cluster IDs for which they are markers. Has to have "gene" and "cluster" as column names
#' @param logFC LogFC threshold for extracting markers from reference clusters using MAST 
#' @param use_reference_for_weights Logical to use either reference or target data for sorting the signature genes
#' @param likelihood_threshold Minimum required likelihood of gene signature score for a cell to be assigned to the respective reference cluster
#' @return list of labels for target cells
#' @return table of marker genes per cluster
#' @export
scid_multiclass <- function(target_gem=NULL, reference_gem=NULL, reference_clusters=NULL, 
                             logFC=0.5, use_reference_for_weights=FALSE, likelihood_threshold=0.99, 
                             markers=NULL) {
  #----------------------------------------------------------------------------------------------------
  # Data preprocessing
  
  # Reference
  # Check if one of reference data (gem and labels) or markers have been given
  if (is.null(reference_gem) && is.null(reference_clusters) && is.null(markers)) {
    message("Please provide either clustered reference data or list of markers for each reference cluster")
    return()
  } else if (is.null(markers)) {
    # Check all reference cells have labels
    if (is.null(reference_clusters)) {
      message("None of the reference cells has a cluster ID. Please check the reference_clusters list provided.")
      return()
    } else {
      common_cells <- intersect(names(reference_clusters), colnames(reference_gem))
      if (length(common_cells) == 0) {
        message("None of the reference cells has a cluster ID. Please check the reference_clusters list provided.")
        return()
      } else {
        reference_gem <- reference_gem[, common_cells]
        rownames(reference_gem) <- toupper(rownames(reference_gem))
        reference_clusters <- reference_clusters[common_cells]
      }
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
  
  # ----------------------------------------------------------------------------------------------------
  # Stage 1: Find signature genes from reference data
  if (is.null(markers)) {
    # ----------------------------------------------------------------------------------------------------
    # Stage 1: Find signature genes from reference data
    message("Stage 1: extract signatures genes from reference clusters")
    markers <- find_markers(reference_gem, reference_clusters, logFC)
    # Filter out signature genes that are not present in the target data
    markers <- markers[which(markers$gene %in% rownames(target_gem)), ]
    celltypes <- unique(markers$cluster)
  } else {
    markers <- markers[which(markers$gene %in% rownames(target_gem)), ]
    celltypes <- unique(markers$cluster)
  }

  # ----------------------------------------------------------------------------------------------------
  # Min-max normalization of target gem
  target_gem_norm <- t(apply(target_gem[markers$gene, ], 1, function(x) normalize_gem(x)))
  na.values <- which(apply(target_gem_norm, 1, function(x){any(is.na(x))}))
  if (length(na.values) > 0) {
    target_gem_norm <- target_gem_norm[-na.values, ]
  }
  
  # ----------------------------------------------------------------------------------------------------
  # Stage 2: Weight signature genes
  message("Stage 2: Estimate weights of signature genes")
  
  weights <- list()
  
  if (use_reference_for_weights) {
    if (is.null(reference_gem)) {
      message("Reference data not available. Please provide clustered reference data or set use_reference_for_weights=FALSE.")
      return()
    }
    reference_gem_norm <- t(apply(reference_gem[markers$gene, ], 1, function(x) normalize_gem(x)))
    na.values <- which(apply(reference_gem_norm, 1, function(x){any(is.na(x))}))
    if (length(na.values) > 0) {
      reference_gem_norm <- reference_gem_norm[-na.values, ]
    }
    for (i in 1:length(celltypes)) {
      svMisc::progress(i*100/length(celltypes))
      Sys.sleep(1 / length(celltypes))
      signature_genes <- markers$gene[which(markers$cluster == celltypes[i])]
      IN <- names(which(reference_clusters == celltypes[i]))
      OUT <- setdiff(colnames(reference_gem), IN)
      weights[[as.character(celltypes[i])]] <- scID_weight(reference_gem_norm[signature_genes, ], IN, OUT)
      if (i==length(celltypes)) cat("Done!")
    }
  } else {
    for (i in 1:length(celltypes)) {
      svMisc::progress(i*100/length(celltypes))
      Sys.sleep(1 / length(celltypes))
      signature_genes <- markers$gene[which(markers$cluster == celltypes[i])]
      putative_groups <- choose_unsupervised(target_gem[markers$gene, ], signature_genes)
      weights[[as.character(celltypes[i])]] <- scID_weight(target_gem_norm[signature_genes, ], putative_groups$in_pop, putative_groups$out_pop)
      if (i==length(celltypes)) cat("Done!")
    }
  }
  
  #----------------------------------------------------------------------------------------------------
  # Stage 3
  # Find scores and putative matches
  message("Stage 3.1-2: Calculate scores and find matching cells")
  
  scores <- data.frame(matrix(NA, length(celltypes), ncol(target_gem)), row.names = celltypes)
  colnames(scores) <- colnames(target_gem)
  
  for (i in 1:length(celltypes)) {
    celltype <- celltypes[i]
    # Remove progress because it is overlapping with mixtools messages
    # svMisc::progress(i*100/length(celltypes))
    # Sys.sleep(1/length(celltypes))
    signature <- names(weights[[as.character(celltypes[i])]])
    weighted_gem <- weights[[as.character(celltypes[i])]] * target_gem_norm[signature, ]
    score <- colSums(weighted_gem)/sum(weights[[as.character(celltypes[i])]])
    matches <- final_populations(score, likelihood_threshold)
    scores[as.character(celltype), matches] <- scale(score[matches])
    if (i==length(celltypes)) cat("Done!")
  }
  
  scID_labels <- c()
  for (i in 1:ncol(scores)) {
    cell <- colnames(scores)[i]
    svMisc::progress(i*100/ncol(scores))
    Sys.sleep(1/ncol(scores))
    if (all(is.na(scores[, cell]))) {
      matching_type <- "unassigned"
    } else {
      matching_type <- rownames(scores)[which(scores[, cell] == max(scores[, cell], na.rm = T))]
    }
    scID_labels <- c(scID_labels, matching_type)
  }
  names(scID_labels) <- colnames(scores)
  
  # return result
  return(list(labels=scID_labels, markers=markers))
}

#' Main function to get Gene expression matrix and signature genes and return matches and scores
#' @param target_gem Data frame of gene expression (rows) per cell (columns) in target data
#' @param gene_signature List of genes symbols enriched in the population of interest
#' @param likelihood_threshold Minimum required likelihood of gene signature score for a cell to be assigned to the respective reference cluster
#' @return list of target cells matching to the population of interest
#' @return list of matching scores for all target cells
#' @export
scid_singleclass <- function(target_gem=NULL, gene_signature=NULL, likelihood_threshold=0.99) {
  #----------------------------------------------------------------------------------------------------
  # Data preprocessing
  
  rownames(target_gem) <- toupper(rownames(target_gem))
  gene_signature <- intersect(toupper(gene_signature), rownames(target_gem))
  
  # ----------------------------------------------------------------------------------------------------
  # Min-max normalization of target gem
  target_gem_norm <- t(apply(target_gem[gene_signature, ], 1, function(x) normalize_gem(x)))
  na.values <- which(apply(target_gem_norm, 1, function(x){any(is.na(x))}))
  if (length(na.values) > 0) {
    target_gem_norm <- target_gem_norm[-na.values, ]
  }
  
  # ----------------------------------------------------------------------------------------------------
  # Stage 2: Weight signature genes
  message("Stage 2: Estimate weights of signature genes")
  
  putative_groups <- choose_unsupervised(target_gem, gene_signature)
  weights <- scID_weight(target_gem_norm[gene_signature, ], putative_groups$in_pop, putative_groups$out_pop)
  
  #----------------------------------------------------------------------------------------------------
  # Stage 3
  # Find scores and putative matches
  message("Stage 3.1-2: Calculate scores and find matching cells")
  
  weighted_gem <- weights * target_gem_norm[names(weights), ]
  score <- colSums(weighted_gem)/sum(weights)
  matches <- final_populations(score, likelihood_threshold)

  return(list(scores=score, matches=matches))
}
