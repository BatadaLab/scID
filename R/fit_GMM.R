#' @export
# final_populations <- function(score, likelihood_threshold) {
#   results <-tryCatch(mixtools::normalmixEM(score, k=2, fast = TRUE, verb = FALSE),error = function(e) {return(NULL)})
#   if (is.null(results)) {
#     return(c())
#   } else {
#     mixmdl <- mixtools::normalmixEM(score, k=2, maxit = 10000, fast = TRUE)
#     # Ask user to choose allowed contamination percentage
#     matches <- intersect(names(score)[which(mixmdl$posterior[, which(mixmdl$mu == max(mixmdl$mu))] >= likelihood_threshold)], 
#                     names(score)[which(score >= min(mixmdl$mu))])
#     
#     return(matches)
#   }
# }

final_populations <- function(score) {
  
  fit <- mclust::densityMclust(score)
  
  # Calculate average scID score per group of cells
  avgScore <- rep(NA, length(unique(fit$classification)))
  names(avgScore) <- unique(fit$classification)
  for (ID in names(avgScore)) {
    avgScore[ID] <- mean(score[names(which(fit$classification == ID))])
  }
  
  matches <- names(fit$classification)[which(fit$classification == names(which(avgScore == max(avgScore))))]
  
}



