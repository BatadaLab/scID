#' Function to row-normalize a numerical vector by the 99th percentile 
#' @param x Numerical vector
#' @return Vector with normalized values
#' @export
normalize_gene <- function(x) {
  
  x <- x+1 
  # Added na.rm=TRUE -- see if changes anything
  pctl <- quantile(x, 0.99, na.rm = TRUE)
  x_norm <- x / pctl
  x_norm[which(x_norm > 1)] <- 1
  
  return(x_norm)
}

#' Function to convert gene counts CPM-normalized
#' @param gem Data frame of 
#' @return CPM-normalized gene expression matrix 
#' @export
counts_to_cpm <- function (counts_gem) {
  return(t(t(counts_gem)/colSums(counts_gem))) * 10^6
}
