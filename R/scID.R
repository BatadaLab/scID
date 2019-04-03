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
                             logFC = 0.5, estimate_weights_from_target = FALSE, weights = NULL, only_pos=FALSE){
  
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
    if (estimate_weights_from_target) {
      rm(reference_gem, reference_clusters)
    }
  } else {
    markers <- markers[which(markers$gene %in% rownames(target_gem)), ]
    celltypes <- unique(markers$cluster)
  }

  # ----------------------------------------------------------------------------------------------------
  # Min-max normalization of target gem
  target_gem_norm <- t(apply(target_gem[unique(markers$gene), ], 1, function(x) normalize_gem(x)))
  na.values <- which(apply(target_gem_norm, 1, function(x){any(is.na(x))}))
  if (length(na.values > 0)) {
    target_gem_norm <- target_gem_norm[-na.values, ]
  }
  
  # ----------------------------------------------------------------------------------------------------
  # Stage 2: Weight signature genes
  if (is.null(weights)) {
    if (estimate_weights_from_target) {
      print("Stage 2: Estimate weights of signature genes from target")
      weights <- list()
      # training_labels <- data.frame(matrix(NA, nrow = ncol(target_gem), ncol = length(celltypes)), row.names = colnames(target_gem))
      # colnames(training_labels) <- celltypes
      
      for (i in 1:length(celltypes)) {
        celltype_markers <- markers[which(markers$cluster == celltypes[i]), ]
        positive_markers <- celltype_markers$gene[which(celltype_markers$avg_logFC > 0)]
        negative_markers <- celltype_markers$gene[which(celltype_markers$avg_logFC < 0)]
        
        putative_groups <- choose_unsupervised(target_gem, positive_markers, negative_markers)
        # training_labels[putative_groups$in_pop, i] <- "IN"
        # training_labels[putative_groups$out_pop, i] <- "OUT"
        signature_genes <- c(positive_markers, negative_markers)
        gene.weights <- scID_weight(target_gem_norm[signature_genes, ], putative_groups$in_pop, putative_groups$out_pop)
        weights[[as.character(celltypes[i])]] <- gene.weights
        svMisc::progress(i*100/length(celltypes), max.value = 100, char = "-", progress.bar = T)
        Sys.sleep(0.01)
        if (i==length(celltypes)) cat("Done!")
      }
      names(weights) <- celltypes
    } else {
      if (!is.null(reference_gem) && !is.null(reference_clusters)) {
        print("Stage 2: Estimate weights of signature genes from reference")
        # training_labels <- NULL
        weights <- list()
        # Normalize reference gem
        ref_gem_norm <-  t(apply(reference_gem[unique(markers$gene), ], 1, function(x) normalize_gem(x)))
        na.values <- which(apply(ref_gem_norm, 1, function(x){any(is.na(x))}))
        if (length(na.values > 0)) {
          ref_gem_norm <- ref_gem_norm[-na.values, ]
        }
        for (i in 1:length(celltypes)) {
          celltype_markers <- markers[which(markers$cluster == celltypes[i]), ]
          positive_markers <- celltype_markers$gene[which(celltype_markers$avg_logFC > 0)]
          negative_markers <- celltype_markers$gene[which(celltype_markers$avg_logFC < 0)]
          signature_genes <- c(positive_markers, negative_markers)
          in_pop <- names(which(reference_clusters == celltypes[i]))
          out_pop <- names(which(reference_clusters != celltypes[i]))
          gene.weights <- scID_weight(ref_gem_norm[signature_genes, ], in_pop, out_pop)
          weights[[as.character(celltypes[i])]] <- gene.weights
          svMisc::progress(i*100/length(celltypes), max.value = 100, char = "-", progress.bar = T)
          Sys.sleep(0.01)
        }
        # Won't need reference data any more
        rm(reference_gem, reference_clusters, ref_gem_norm)
      } else {
        print("Please provide reference data in order to calculate weights. Alternative select to estimate weights from target data or provide precompted gene weights.")
        return()
      }
    }
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
    # scores[as.character(celltype), ] <- score
    matches <- final_populations(score) 
    scores[as.character(celltype), matches] <- scale(score[matches])
    if (i==length(celltypes)) cat("Done!")
  }
  # return(list(scores=scores, markers=markers, weights=weights))
  
  # Resolve multiclass assignments
  print ("Stage 3.3: Resolve multiclass assignments")
  labels <- apply(scores, 2, function(x) {ifelse(all(is.na(x)), "unassigned", rownames(scores)[which(x == max(x, na.rm = T))])})

  # return result
  return(list(scores=scores, markers=markers, estimated_weights=weights, labels=labels))

}
