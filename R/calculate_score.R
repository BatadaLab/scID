#' @export
final_populations <- function(score, contamination) {
  results <- capture.output(mixtools::normalmixEM(score))
  if (results[1] %in% c("WARNING! NOT CONVERGENT! ", "Error in normalmixEM(matching_score) : Too many tries! ")) {
    matches <- c()
  } else {
    mixmdl <- mixtools::normalmixEM(score, k=2, maxit = 10000, fast = TRUE)
    # Ask user to choose allowed contamination percentage
    IN <- intersect(names(score)[which(mixmdl$posterior[, which(mixmdl$mu == max(mixmdl$mu))] >= 1-contamination)], 
                    names(score)[which(score >= min(mixmdl$mu))])
    OUT <- intersect(names(score)[which(mixmdl$posterior[, which(mixmdl$mu == min(mixmdl$mu))] > 0.5)], 
                     names(score)[which(score <= min(mixmdl$mu))])
    
    return(list(IN=IN, OUT=OUT))
  }
}

#' @export
adjust_score <- function(scData, matching_score, hk_genes, dtstID=NULL) {
  # Fit linear regression between score for signature and average expression of housekeeping genes
  if (!is.null(dtstID)) {
    classes <- unique(dtstID)
    adjusted_score <- c()
    for (i in 1:length(classes)) {
      gem <- scData[, names(which(dtstID == classes[i]))]
      avg_hk_expression <- colMeans(gem[intersect(rownames(gem), hk_genes), ])
      data <- data.frame(matching_score=matching_score, avg_hk_expression=avg_hk_expression[names(matching_score)])
      fit <- lm(matching_score ~ avg_hk_expression, data = data)
      adjusted_score <- c(adjusted_score, residuals(fit))
    }
    # Additionally, adjust for different number of genes per cell between experiments
    genes_per_cell <- apply(scData, 2, function(x) sum(x>0))
    data <- data.frame(matching_score=adjusted_score, genes_per_cell=genes_per_cell)
    fit <- lm(matching_score ~ genes_per_cell, data = data)
    adjusted_score <- residuals(fit)
  } else {
    avg_hk_expression <- colMeans(scData[intersect(rownames(scData), hk_genes), ])
    data <- data.frame(matching_score=matching_score, avg_hk_expression=avg_hk_expression[names(matching_score)])
    fit <- lm(matching_score ~ avg_hk_expression, data = data)
    adjusted_score <- residuals(fit)
  }
  return(adjusted_score)
}




