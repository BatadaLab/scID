#' Function to row-normalize a numerical vector by the 99th percentile 
#' 
#' @param x Numerical vector
#' 
#' @return Vector with normalized values
#' @export
normalize_gene <- function(x) {
  
  x <- x + 1 
  pctl <- quantile(x, 0.99, na.rm = TRUE)
  x_norm <- pmin(x / pctl, 1)
  
  x_norm
}

#' Function to convert gene counts CPM-normalized
#' @param gem Data frame of 
#' @return CPM-normalized gene expression matrix 
#' @export
counts_to_cpm <- function (counts_gem) {
  
  # Discard genes that are zero across all cells
  IDX_ZEROS <- apply(counts_gem, 1, function(row) all(row==0))
  counts_gem <- counts_gem[!IDX_ZEROS, ]
  print(sprintf("After discarding all zero rows, left with %s rows.", dim(counts_gem)[1]))
  
  counts_scater <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(counts_gem)))
  gem_cpm = scater::calculateCPM(counts_scater)
  
  gem_cpm
}
