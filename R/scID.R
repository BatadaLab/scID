#' Main function to get Gene expression matrix and signature genes and return matches and scores
#' @param target_gem Data frame of gene expression (rows) per cell (columns) in target data
#' @param reference_gem Data frame of gene expression (rows) per cell (columns) in reference data
#' @param reference_clusters Named list of cluster IDs of reference cells
#' @param markers Data frame of cluster specific genes that will be used instead of the reference data
#' @param logFC LogFC threshold for extracting markers from reference clusters
#' @param use_reference_for_weights Logical to use either reference or target data for sorting the signature genes
#' @param only_pos Logical to include negative markers in the cluster specific gene sets
#' @param normalize_reference Logical to select if reference data need to be normalized (when raw counts have been provided)
#' @return list of cluster IDs for target cells
#' @return scID scores for target cells 
#' @return markers used 
#' @return estimated gene weights
#' @export
scid_multiclass <- function(target_gem = NULL, reference_gem = NULL, reference_clusters = NULL, markers = NULL,
                            logFC = 0.5, estimate_weights_from_target = FALSE, weights = NULL, only_pos=FALSE, normalize_reference=TRUE){
  
  #----------------------------------------------------------------------------------------------------
  # Data pre-processing
  if (is.null(reference_gem) && is.null(reference_clusters) && is.null(markers)) {
    stop("Please provide either clustered reference data or list of markers for each reference cluster")
  } else if (is.null(markers)) {
    # Check all reference cells have a cluster ID
    common_cells <- intersect(names(reference_clusters), colnames(reference_gem))
    if (length(common_cells) == 0) {
      stop("None  of the reference cells has a cluster ID. Please check the reference_clusters list provided.")
    } else {
      reference_gem <- reference_gem[, common_cells]
      rownames(reference_gem) <- make.names(toupper(rownames(reference_gem)), unique=TRUE)
        
      # Remove genes that are zero across all cells
      reference_gem <- reference_gem[which(rowSums(reference_gem) != 0), ]
      reference_clusters <- reference_clusters[common_cells]
    }
  } else {
    # Check markers have gene and cluster columns
    if (length(intersect(c("gene", "cluster"), colnames(markers))) !=2 ) {
      stop("Please provide a data frame of markers with gene and cluster in columns")
    }
    if (!"gene" %in% colnames(markers)) {
      stop("Please provide a data frame of markers with gene and cluster in columns")
    } else {
      markers$gene <- toupper(markers$gene)
    }
  }

  # Target
  rownames(target_gem) <- make.names(toupper(rownames(target_gem)), unique=TRUE)
  # Remove genes that are zero across all cells
  target_gem <- target_gem[which(rowSums(target_gem) != 0), ]

  # ----------------------------------------------------------------------------------------------------
  # Stage 1: Find signature genes from reference data
  if (is.null(markers)) {
    # ----------------------------------------------------------------------------------------------------
    # Stage 1: Find signature genes from reference data
    message("Stage 1: extract signatures genes from reference clusters")
    markers <- find_markers(reference_gem, reference_clusters, logFC, only.pos=only_pos, normalize_reference=normalize_reference)
    # Filter out signature genes that are not present in the target data
    markers <- markers[which(markers$gene %in% rownames(target_gem)), ]
    celltypes <- unique(markers$cluster)
    if (estimate_weights_from_target) {
      rm(reference_gem, reference_clusters)
    }
  } else {
    markers <- markers[which(markers$gene %in% rownames(target_gem)), ]
    celltypes <- unique(markers$cluster)
  }

  # ----------------------------------------------------------------------------------------------------
  # Min-max normalization of target gem
  target_gem_norm <- t(apply(target_gem[unique(markers$gene), ], 1, function(x) normalize_gene(x)))
  target_gem_norm <- target_gem_norm[complete.cases(target_gem_norm), ]
  
  # ----------------------------------------------------------------------------------------------------
  # Stage 2: Weight signature genes
  if (is.null(weights)) {
    if (estimate_weights_from_target) {
      message("Stage 2: Estimate weights of signature genes from target")
      weights <- list()
      for (i in 1:length(celltypes)) {
        svMisc::progress(i*100/length(celltypes), max.value = 100, char = "-", progress.bar = T)
        Sys.sleep(0.01)
        celltype_markers <- markers[which(markers$cluster == celltypes[i]), ]
        positive_markers <- celltype_markers$gene[which(celltype_markers$avg_logFC > 0)]
        negative_markers <- celltype_markers$gene[which(celltype_markers$avg_logFC < 0)]
        training_groups <- choose_training_set(target_gem, positive_markers, negative_markers)
        signature_genes <- c(positive_markers, negative_markers)
        gene.weights <- scID_weight(target_gem_norm[signature_genes, , drop=FALSE], training_groups$in_pop, training_groups$out_pop)
        weights[[as.character(celltypes[i])]] <- gene.weights
        if (i==length(celltypes)) cat("Done!")
      }
      names(weights) <- celltypes
    } else {
      if (!is.null(reference_gem) && !is.null(reference_clusters)) {
        message("Stage 2: Estimate weights of signature genes from reference")
        weights <- list()
        # Normalize reference gem
        ref_gem_norm <-  t(apply(reference_gem[unique(markers$gene), ], 1, function(x) normalize_gene(x)))
        ref_gem_norm <- ref_gem_norm[complete.cases(ref_gem_norm), ]
        for (i in 1:length(celltypes)) {
          svMisc::progress(i*100/length(celltypes))#, max.value = 100, char = "-", progress.bar = T)
          Sys.sleep(0.01)
          signature_genes <- markers[which(markers$cluster == celltypes[i]), "gene"]
          true_cells <- names(reference_clusters)[which(reference_clusters == as.character(celltypes[i]))]
          false_cells <- setdiff(names(reference_clusters), true_cells)
          gene.weights <- scID_weight(gem = ref_gem_norm[signature_genes, ,drop=FALSE], true_cells, false_cells)
  
          weights[[as.character(celltypes[i])]] <- gene.weights
          if (i==length(celltypes)) cat("Done!")
        }
        # Won't need reference data any more, remove for efficiency
        rm(reference_gem, reference_clusters, ref_gem_norm)
      } else {
        stop("Please provide reference data in order to calculate weights, choose to estimate weights from target data, or provide precompted gene weights.")
      }
    }
  }

  #----------------------------------------------------------------------------------------------------
  # Stage 3
  # Find scores and putative matches
  message("Stage 3.1-2: Calculate scores and find matching cells")
  
  scores <- data.frame(matrix(NA, length(celltypes), ncol(target_gem)), row.names = celltypes)
  colnames(scores) <- colnames(target_gem)
  
  full_scores <- data.frame(matrix(NA, length(celltypes), ncol(target_gem)), row.names = celltypes)
  colnames(full_scores) <- colnames(target_gem)
  
  for (i in 1:length(celltypes)) {
    celltype <- as.character(celltypes[i])
    signature <- intersect(names(weights[[celltype]]), rownames(target_gem_norm))
    weighted_gem <- weights[[celltype]][signature] * target_gem_norm[signature, ,drop=FALSE]
    
    score <- colSums(weighted_gem)/sqrt(sum(weights[[celltype]]^2))
    matches <- final_populations(score) 
    scores[as.character(celltype), matches] <- scale(score[matches])
    full_scores[as.character(celltype), ] <- score
    
    if (i==length(celltypes)) cat("Done!")
  }
  
  # Resolve multiclass assignments
  message ("Stage 3.3: Resolve multiclass assignments")
  labels <- apply(scores, 2, function(x) {ifelse(all(is.na(x)), "unassigned", rownames(scores)[which(x == max(x, na.rm = T))])})

  # return result
  return(list(markers=markers, estimated_weights=weights, labels=labels, scores=full_scores))

}

