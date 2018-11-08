#' @export
final_populations <- function(score, contamination, lambda_est=NULL, mean_est=NULL, sd_est=NULL) {
  results <-tryCatch(mixtools::normalmixEM(score, k=2, fast = TRUE, verb = FALSE),error = function(e) {return(NULL)})
  if (is.null(results)) {
    return(c())
  } else {
    if (is.null(lambda_est) && is.null(mean_est) && is.null(sd_est)) {
      mixmdl <- mixtools::normalmixEM(score, k=2, maxit = 10000, fast = TRUE)
    } else {
      mixmdl <- mixtools::normalmixEM(score, k=2, maxit = 10000, fast = TRUE, lambda = lambda_est, mu = mean_est, sigma = sd_est)
    }
    # Ask user to choose allowed contamination percentage
    IN <- intersect(names(score)[which(mixmdl$posterior[, which(mixmdl$mu == max(mixmdl$mu))] >= 1-contamination)], 
                    names(score)[which(score >= min(mixmdl$mu))])
    OUT <- intersect(names(score)[which(mixmdl$posterior[, which(mixmdl$mu == min(mixmdl$mu))] > 0.5)], 
                     names(score)[which(score <= min(mixmdl$mu))])
    
    return(IN)
  }
}





