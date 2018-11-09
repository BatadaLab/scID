#' @export
normalize_gem <- function(x) {
  
  x <- x+1 
  # Added na.rm=TRUE -- see if changes anything
  pctl <- quantile(x, 0.99, na.rm = TRUE)
  x_norm <- x / pctl
  x_norm[which(x_norm > 1)] <- 1
  
  return(x_norm)
}

