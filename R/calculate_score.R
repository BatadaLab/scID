#' @export

final_populations <- function(score, contamination) {
  results <- capture.output(mixtools::normalmixEM(score))
  if (results[1] %in% c("WARNING! NOT CONVERGENT! ", "Error in normalmixEM(matching_score) : Too many tries! ")) {
    matches <- c()
  } else {
    mixmdl <- mixtools::normalmixEM(score, k=2, maxit = 10000, fast = TRUE)
    #plot(mixmdl,which=2)
    # select threshold based on purity
    if (mixmdl$mu[1] > mixmdl$mu[2]) {
      j <- 1
      i <-2
    } else {
      j <- 2
      i <- 1
    }
    curve <- mixmdl$lambda[j] * dnorm(score, m=mixmdl$mu[j], sd=mixmdl$sigma[j]) / (mixmdl$lambda[j] * dnorm(score, m=mixmdl$mu[j], sd=mixmdl$sigma[j]) + mixmdl$lambda[i] * dnorm(score, m=mixmdl$mu[i], sd=mixmdl$sigma[i]))
    #plot(matching_score, curve)
    
    # Ask user to choose allowed contamination percentage
    IN <- names(which(curve >= (1-contamination)))
    OUT <- names(which(curve <= quantile(curve, 0.1)))
    
    return(list(IN=IN, OUT=OUT))
  }
}

adjust_score <- function(scData, matching_score, hk_genes) {
  
  # Fit linear regression between score for signature and average expression of housekeeping genes
  avg_hk_expression <- colMeans(scData[intersect(rownames(scData), hk_genes), ])
  
  data <- data.frame(matching_score=matching_score, avg_hk_expression=avg_hk_expression[names(matching_score)])
  fit <- lm(matching_score ~ avg_hk_expression, data = data)
  adjusted_score <- residuals(fit)
  
  return(adjusted_score)
}




