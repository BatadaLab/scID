#' Function to identify matching population from scID score by fitting a Gaussian finite mixture model from Mclust
#' @param score List of scID scores for target cells
#' @return List of matching cells
#' @export

final_populations <- function(score) {
  
  fit <- mclust::densityMclust(score)
  
  # Calculate average scID score per group of cells
  avgScore <- rep(NA, length(unique(fit$classification)))
  names(avgScore) <- unique(fit$classification)
  for (ID in names(avgScore)) {
    avgScore[ID] <- mean(score[names(which(fit$classification == ID))])
  }
  
  matches <- names(fit$classification)[which(fit$classification == names(which(avgScore == max(avgScore))))]
  
  return(matches)
}
