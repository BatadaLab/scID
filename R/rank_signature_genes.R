#' Function to choose initial IN and OUT populations with linear regression
#' @param gem Data frame of gene expression of genes (rows) in cells (columns)
#' @param signature_genes list of gene names (should be same format as rownames of gem)
#' @return a lists of confident IN and OUT cells
#' @export
choose_unsupervised <- function(gem, signature_genes) {
  
  # Bin values to 0 and 1 for present (expressed) and absent genes
  binned_gem <- apply(gem, 1, function(x) ifelse(x>quantile(x[which(x>0)], 0.25, na.rm = TRUE), 1, 0))
  
  # Find total number of expressed genes per cell (n_e)
  n_e <- apply(binned_gem, 1, function(x) length(which(x==1)))
  # Find total number of expressed marker genes per cell (n_me)
  n_me <- apply(binned_gem[, which(colnames(binned_gem) %in% signature_genes)], 1, function(x) length(which(x==1)))
  # Find total number of marker genes (n_m)
  n_m <- length(signature_genes)
  
  data <- data.frame(recall = n_me/n_m, precision = n_me/n_e)
  
  na.values <- which(apply(data, 1, function(x){any(is.na(x))}))
  
  if (length(na.values) > 0) {
    data <- data[-na.values, ]
  }
  
  clustered_data <- kmeans(scale(data), 2, iter.max = 1000)
  
  clustered_data$centers
  IN_id <- intersect(which(clustered_data$centers[, "recall"] == max(clustered_data$centers[, "recall"])), 
                     which(clustered_data$centers[, "precision"] == max(clustered_data$centers[, "precision"])))
  
  return(list(in_pop=names(which(clustered_data$cluster == IN_id)), out_pop=names(which(clustered_data$cluster != IN_id))))
}

#' Main function for estimation of gene ranks
#' @param gem Data frame of signature genes in cells
#' @param true_cells list of cell names that are confidently IN-population
#' should match with column names of gem
#' @param false_cells list of cell names that are confidently OUT-population
#' should match with column names of gem
#' @return list of weights for each signature gene
#' @export
scID_weight <- function(gem, true_cells, false_cells) {
  
  weights <- rep(NA, nrow(gem))
  names(weights) <- rownames(gem)
  for (gene in rownames(gem)) {
    numerator <- mean(as.numeric(gem[gene, true_cells])) - mean(as.numeric(gem[gene, false_cells]))
    denominator <- sd(as.numeric(gem[gene, false_cells]))^2 + sd(as.numeric(gem[gene, true_cells]))^2
    weights[gene] <- max(numerator/denominator, 0)
  }
  
  weights[which(is.na(weights))] <- 0
  
  return(weights/max(weights))
}

