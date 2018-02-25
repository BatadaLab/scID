#' @export

#' Function to choose initial IN and OUT populations with linear regression
#' @param gem Data frame of gene expression of genes (rows) in cells (columns)
#' @param signature_genes list of gene names (should be same format as rownames of gem)
#' @return a lists of confident IN and OUT cells
choose_unsupervised <- function(gem, signature_genes) {
  cellPct <- c()
  binned_gem <- apply(gem, 1, function(x) ifelse(x>quantile(x[which(x>0)], 0.25, na.rm = TRUE), 1, 0))
  
  all_gene_counts <- apply(binned_gem, 1, function(x) length(which(x==1)))
  signature_gene_counts <- apply(binned_gem[, which(colnames(binned_gem) %in% signature_genes)], 1, function(x) length(which(x==1)))
  
  data <- data.frame(signature_gene_counts=signature_gene_counts, all_gene_counts=all_gene_counts)
  fit <- lm(signature_gene_counts ~ all_gene_counts, data = data)
  intersept <- fit$coefficients[1]
  slope <- fit$coefficients[2]
  predicted_scores <- c()
  for (i in 1:nrow(data)) {
    predicted_scores[i] <- intersept + (slope * data[i,"all_gene_counts"])
  }
  s_r <- sqrt(sum((data[, "signature_gene_counts"] - predicted_scores)^2) / (nrow(data) - 2))
  
  positive_thr <- predicted_scores + 2 * s_r
  names(positive_thr) <- rownames(data)
  true_differences <- data[, "signature_gene_counts"] - positive_thr
  true_cells <- names(which(true_differences > 0))
  
  negative_thr <- predicted_scores - s_r
  names(negative_thr) <- rownames(data)
  false_differences <- data[, "signature_gene_counts"] - negative_thr
  false_cells <- names(which(false_differences < 0))
  
  if (length(false_cells) == 0) {
    remaining_cells = setdiff(colnames(gem), true_cells)
    false_cells = remaining_cells[sample(length(remaining_cells), min(2*length(true_cells), 0.5*length(remaining_cells)))]
  }
  
  return(list(in_pop=true_cells, out_pop=false_cells))
}

#' Function to select initial IN population from necessary positive markers
#' @param gem Data frame of expression of positive/necessary markers (rows) in cells of 
#' the dataset (columns)
#' @return lists of confident IN (cells that express all the positive markers)
#' and confident OUT cells (cells that do not express any of the positive markers)
choose_with_positive_markers <- function(gem) {
  expressing_cells <- apply(gem, 1, function(x) names(which(x>quantile(x[which(x > 0)], 0.2, na.rm = TRUE))))
  nonexpressing_cells <- apply(gem, 1, function(x) names(which(x==0)))
  if (nrow(gem) == 1) {
    true_cells <- as.vector(expressing_cells)
    false_cells <- as.vector(nonexpressing_cells)
  } else {
    true_cells <- Reduce(intersect, expressing_cells)
    false_cells <- Reduce(intersect, nonexpressing_cells)
  }
  if (length(false_cells) == 0) {
    return(TRUE)
  } else if (length(true_cells) == 0) {
    return(FALSE)
  } else {
    return(list(in_pop=true_cells, out_pop=false_cells))
  }
}

#' Function to select initial OUT population from negative markers
#' @param gem Data frame of expression of negative markers (rows) in cells of the dataset (columns)
#' @return list of confident OUT population
choose_with_negative_markers <- function(gem) {
  expressing_cells <- apply(gem, 1, function(x) names(which(x>quantile(x[which(x > 0)], 0.7, na.rm = TRUE))))
  if (nrow(gem) == 1) {
    false_cells <- as.vector(expressing_cells)
  } else {
    false_cells <- Reduce(union, expressing_cells)
  }
  return(false_cells)
}

#' Function to calculate Signal-to-noise Ratio between IN and OUT populations
#' @param gem Data frame of signature genes in cells
#' @param true_cells list of cell names that are confidently IN-population
#' should match with column names of gem
#' @param false_cells list of cell names that are confidently OUT-population
#' should match with column names of gem
#' @return list of SNR value for each signature gene
calculate_snr <- function(gem, true_cells, false_cells) {
  
  snr <- rep(NA, nrow(gem))
  names(snr) <- rownames(gem)
  for (gene in rownames(gem)) {
    signal <- mean(as.numeric(gem[gene, true_cells]))
    noise <- sd(as.numeric(gem[gene, false_cells]))
    snr[gene] <- signal/(noise+1) 
  }
  
  return(snr)
}

#' Main function for estimation of gene ranks
#' @param gem Data frame of gene expression (rows) in cells (columns)
#' @param signature List of signature gene names (should match with rownames of gem)
#' @param positive_markers (optional) List of genes that should be present in all confident IN-cells
#' @param negative_markers (optional) List of genes that should not be present in any IN-cell
#' @return List of weights per gene
weight_signature <- function(gem, signature, positive_markers=NULL, negative_markers=NULL) {
  
  if (is.null(positive_markers)) {
    res <- choose_unsupervised(gem, signature)
  } else {
    res <- choose_with_positive_markers(na.omit(gem[positive_markers, ]))
  }
  
  if (typeof(res)=="logical") {
    if (res==TRUE) {
      return("All cells express the positive markers")
    } else if (res==FALSE) {
      return("No cell expresses the positive markers")
    }
  } else {
    true_cells <- res$in_pop
    if (length(true_cells) == 0) {
      genes_per_cell <- apply(na.omit(gem[signature, ]), 2, function(x) sum(x>0))
      true_cells <- names(which(genes_per_cell > quantile(genes_per_cell, 0.77)))
    }
    if (!is.null(negative_markers)) {
      false_cells <- choose_with_negative_markers(na.omit(gem[negative_markers, ]))
    } else {
      false_cells <- res$out_pop
    }
    common_cells <- intersect(true_cells, false_cells)
    if (length(common_cells) > 0) {
      true_cells <- setdiff(true_cells, common_cells)
      false_cells <- setdiff(false_cells, common_cells)
    }
    gem_sig <- na.omit(gem[signature, c(true_cells, false_cells)])
    gene_correlations <- cor(t(gem_sig), method = "spearman")
    nas <- which(rowSums(is.na(gene_correlations))==ncol(gene_correlations)-1)
    if (length(nas) > 0) {
      gene_correlations <- gene_correlations[-nas, -nas]
    }
    mean_cors <- rowMeans(gene_correlations)
    keep <- names(which(mean_cors > 0))
    snr <- calculate_snr(gem[keep, ], true_cells, false_cells)
    weights <- snr/max(snr)
    return(weights) 
  }
}