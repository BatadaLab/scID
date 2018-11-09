#' Main function to get Gene expression matrix and signature genes and return matches and scores
#' @param signature_file (optional) filename of the signature genes
#' @param signature_genes (optional) a list of the signature genes
#' either signature_file or signature_genes must be given as input
#' @param gem_file (optional) filename of the gene expression matrix
#' @param gem (optional) Data frame of gene expression (rows) per cell (columns)
#' either gem_file or gem must be given as input
#' @param contamination Percentage of accepted cells that belong to the common area between 
#' IN and OUT population distributions
#' @param do.imputation Logical to choose not to correct for dropouts (default is TRUE for doing imputation)
#' @param gene.weigts list of weights for the signature genes if known
#' @return list of names of matching cells (IN-population) and 
#' @return matching score of every cell of the dataset
#' @return list of weight per gene
#' @export
scid_match_cells <- function(signature_file=NULL, gem_file=NULL, gem=NULL, signature_genes=NULL, 
                             contamination=0, sort.signature = TRUE, gene.weights=NULL) {
  
  #----------------------------------------------------------------------------------------------------
  # Read data
  if (is.null(gem)) {
    print("Reading data")
    gem <- loadfast(gem_file)
  }
  if (is.null(signature_genes)) {
    print("Reading signature genes")
    signature_genes <- read.table(signature_file, stringsAsFactors = FALSE)$V1
  }
  
  # Make rownames and signature upper case
  rownames(gem) <- toupper(rownames(gem))
  signature_genes <- toupper(signature_genes)
  
  # Keep only genes in scData
  signature_genes <- intersect(signature_genes, rownames(gem))
  
  #----------------------------------------------------------------------------------------------------
  # Find score for signature genes
  # Calculate specificity
  if (length(signature_genes) == 0) {
    print("No gene from the signature is present in the GEM. Please check the given gene list.")
    return(list(matches=NULL, matchingScore=rep(NA, ncol(gem))))
  } else if (length(signature_genes) ==  1) {
    print("Only one gene from the given gene list is present in the GEM. Please check the given list.")
    return(list(matches=NULL, matchingScore=rep(NA, ncol(gem))))
  } else {
    # normalize to 0,1 
    gem_norm <- t(apply(gem[signature_genes, ], 1, function(x) normalize_gem(x)))
    na.values <- which(apply(gem_norm, 1, function(x){any(is.na(x))}))
    if (length(na.values > 0)) {
      gem_norm <- gem_norm[-na.values, ]
    }
    
    if (is.null(weights)) {
      putative_groups <- choose_unsupervised(gem, signature_genes)
      weights <- scID_weight(gem_norm, putative_groups$in_pop, putative_groups$out_pop)
    }
    
    weighted_gem <- weights[rownames(gem_norm)] * gem_norm
    
    score <- colSums(weighted_gem)/sum(weights[rownames(gem_norm)])
    
    if (is.null(weights)) {
      #lambda_est <- c(length(putative_groups$out_pop)/ncol(gem), length(putative_groups$in_pop)/ncol(gem))
      #mean_est <- c(mean(score[putative_groups$out_pop]), mean(score[putative_groups$in_pop]))
      #sd_est <- c(sd(score[putative_groups$out_pop]), sd(score[putative_groups$in_pop]))
      matches <- final_populations(score, contamination)#, lambda_est, mean_est, sd_est)
    } else {
      matches <- final_populations(score, contamination)
    }
    
    return(list(matches=matches, matchingScore=score, weights=weights))
    
  }
}
