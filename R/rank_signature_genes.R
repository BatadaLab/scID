#' Function to choose training IN and OUT populations using precision-recall
#' 
#' @param gem Data frame of gene expression of genes (rows) in cells (columns)
#' @param positive_markers List of gene names expected to be upregulated in IN population
#' @param negative_markers List of gene names expected to be downregulated in IN population
#' 
#' @return Lists of training IN and OUT cells
#' @export
choose_training_set <- function(gem, positive_markers, negative_markers) {
  
  positive_markers <- intersect(positive_markers, rownames(gem))
  negative_markers <- intersect(negative_markers, rownames(gem))
  # Bin values to 0 and 1 for present (expressed) and absent genes
  sink("aux");
  binned_gem <- apply(gem, 1, function(x) biomod2::BinaryTransformation(x, threshold = quantile(x[which(x>0)], 0.25, na.rm = TRUE)))
  sink(NULL);
  
  # Find total number of expressed genes per cell (n_e)
  n_e <- rowSums(binned_gem)
  # Find total number of expressed positive marker genes per cell (n_pme)
  if (length(positive_markers) >= 1) {
    n_pme <- rowSums(binned_gem[, positive_markers, drop = FALSE])
  } else {
    n_pme <- rep(0, nrow(binned_gem))
    names(n_pme) <- rownames(binned_gem)
  }
  # Find total number of expressed negative marker genes per cell (n_nme)
  if (length(negative_markers) >= 1) {
    n_nme <- rowSums(binned_gem[, negative_markers, drop = FALSE])
  } else {
    n_nme <- rep(0, nrow(binned_gem))
    names(n_nme) <- rownames(binned_gem)
  }
  # Find total number of positive marker genes (n_pm)
  n_pm <- max(length(positive_markers), 1)
  # Find total number of negative marker genes (n_nm)
  n_nm <- max(length(negative_markers), 1)
  
  data <- data.frame(
    recall = (n_pme/n_pm) - (n_nme/n_nm), 
    precision = (n_pme-n_nme)/n_e
    )
  rownames(data) <- colnames(gem)
  data <- data[complete.cases(data), ]
  
  library(mclust)
  sink("aux");
  fit <- Mclust(data)
  sink(NULL);
  
  # Get centroids of each cluster
  centroids <- data.frame(matrix(NA, length(unique(fit$classification)), 2), row.names = unique(fit$classification))
  colnames(centroids) <- c("precision", "recall")
  sds <- data.frame(matrix(NA, length(unique(fit$classification)), 2), row.names = unique(fit$classification))
  colnames(sds) <- c("precision", "recall")
  for (ID in rownames(centroids)) {
    centroids[ID, "precision"] <- mean(data[which(fit$classification == ID), "precision"]) 
    sds[ID, "precision"] <- sd(data[which(fit$classification == ID), "precision"]) 
    centroids[ID, "recall"] <- mean(data[which(fit$classification == ID), "recall"]) 
    sds[ID, "recall"] <- sd(data[which(fit$classification == ID), "recall"]) 
  }
  
  IN_candidates <- unique(c(rownames(centroids)[which(centroids$recall == max(centroids$recall))], rownames(centroids)[which(centroids$precision == max(centroids$precision))]))
                    
  E_dist <- apply(centroids, 1, function(x) sqrt((1-x[1])^2 + (1-x[2])^2))
  
  IN_id <- names(E_dist)[which(E_dist == min(E_dist))]
  
  IN_cells <- colnames(gem)[which(fit$classification %in% IN_id)]
  
  # If there are two clusters found as candidate IN remove the one that is farthest from (1,1)
  other_IN <- setdiff(IN_candidates, IN_id)
  if (length(other_IN) == 1) {
    NA_cells <- colnames(gem)[which(fit$classification %in% other_IN)]
  } else {
    NA_cells <- c()
  }
  
  # Get OUT cells removing those that are in the IN radious
  OUT_cells <- setdiff(rownames(data), c(IN_cells, NA_cells))
  
  list(in_pop=IN_cells, out_pop=OUT_cells)
}

#' Main function for estimation of gene ranks
#' 
#' @param gem Data frame of signature genes in cells
#' @param true_cells List of training IN cells
#' @param false_cells List of training OUT cells
#' 
#' @return List of weights for signature genes
#' @export
scID_weight <- function(gem, true_cells, false_cells) {

  weights <- rep(NA, nrow(gem))
  names(weights) <- rownames(gem)
  for (gene in rownames(gem)) {
    numerator <- mean(as.numeric(gem[gene, true_cells])) - mean(as.numeric(gem[gene, false_cells]))
    denominator <- sd(as.numeric(gem[gene, false_cells]))^2 + sd(as.numeric(gem[gene, true_cells]))^2
    
    weights[gene] <- numerator/denominator
  }
  weights[which(is.na(weights))] <- 0

  weights
}
