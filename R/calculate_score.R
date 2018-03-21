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
adjust_score <- function(scData, matching_score, hk_genes) {
  
  # Fit linear regression between score for signature and average expression of housekeeping genes
  avg_hk_expression <- colMeans(scData[intersect(rownames(scData), hk_genes), ])
  
  data <- data.frame(matching_score=matching_score, avg_hk_expression=avg_hk_expression[names(matching_score)])
  fit <- lm(matching_score ~ avg_hk_expression, data = data)
  adjusted_score <- residuals(fit)
  
  return(adjusted_score)
}




