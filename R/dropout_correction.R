#' @export


normalized_exp_pctl <- function(x) {
  
  x <- x+1 
  # Added na.rm=TRUE -- see if changes anything
  pctl <- quantile(x, 0.99, na.rm = TRUE)
  x_norm <- x / pctl
  x_norm[which(x_norm > 1)] <- 1
  
  return(x_norm)
}

#' Initialize vector of paramaters for the mixture model
#' @param x vectors of a gene's expression across all cells
#' @return a named vector with the initial parameters: [lambda, theta, mu, sigma]
initialize_prmt_exp <- function(x) {
  
  threshold <- log2(1.01)
  prmt <- rep(0, 4)
  names(prmt) <- c("lambda", "theta", "mu", "sigma")
  # Initialize lambda
  prmt$lambda <- sum(x == threshold) / length(x)
  # Initialize theta 1/mean(x)
  prmt$theta <- 1/mean(x)
  # Initialize mu and sigma - for this remove zero values (1.01)
  x_nonzero <- x[x > threshold]
  prmt$mu <- mean(x_nonzero)
  prmt$sigma <- sd(x_nonzero)
  if (is.na(prmt$sigma)) { prmt$sigma <- 0 }
  
  return(prmt)
}

#' Update parameters
#' @param d The posterior latent distribution
#' @param x The vector of gene expression across all cells
#' @return the parameters lambda, theta, mu and sigma
update_prmt_exp <- function(d, x) {
  
  # Update lambda
  lambda <- mean(d)
  # Update Normal distribution parameters
  mu <- sum((1-d)*x)/sum(1-d)
  sigma <- sqrt(sum((1-d) * ((x - mu))^2)/sum(1-d))
  # Update exponential distribution parameter
  theta <- length(x) / sum(d*x)
  
  return(list(lambda, theta, mu, sigma))
}

#' Estimate Exponential and Normal distribution parameters with EM algorithm
#' from single-cell data
#' @param gem Gene expression matrix, rows are genes and columns are cells
#' gem should be normalized and log-transformed (log2(gem + 1.01))
#' 1.01 used instead of 1 to avoid zero values
#' @return a table with the estimated dropout probabilities for each gene in each cell of the gem
#' @return a list of estimated parameters for each gene
EM_for_dropout_estimation <- function(gem) {
  
  D <- data.frame(matrix(0, nrow = nrow(gem), ncol = ncol(gem)), row.names = rownames(gem))
  colnames(D) <- colnames(gem)
  
  parameters <- data.frame(matrix(0, nrow = nrow(gem), ncol = 4), row.names = rownames(gem))
  colnames(parameters) <- c("lambda", "theta", "mu", "sigma")
  
  # Estimating mixture models
  threshold <- log2(1.01)
  
  # Find genes that are zero in all cells
  null_genes <- which(abs(rowSums(gem) - threshold * ncol(gem)) < 1e-10)
  for (gene in rownames(gem)) {
    if (gene %in% null_genes) {
      D[gene, ] <- 1 # return all 1 (it is a dropout) -- think about it
      parameters[gene] <- rep(NA, 4)
    } else {
      x <- unlist(gem[gene, ])
      
      # Initialize parameters
      prmt <- initialize_prmt_exp(x)
      eps <- 10
      iter <- 0
      loglik_old <- 0
      
      while(eps > 0.05 & iter <= 500) {
        # Generate Exponential distribution with shape=alpha and rate=beta
        exp_distr = prmt$lambda * dexp(x, rate = prmt$theta, log = FALSE)
        # Genearte Normal distribution with mean=mu and sd=sigma
        norm_distr = (1 - prmt$lambda) * dnorm(x, mean = prmt$mu, sd = prmt$sigma)
        # Equation 4 (d)
        if (identical(exp_distr, numeric(0))) {
          d = 0
        } else {
          d = exp_distr / (exp_distr + norm_distr)
        }
        
        # Update parameters
        new_prmt <- update_prmt_exp(d, x)
        prmt$lambda <- new_prmt[[1]]
        prmt$theta <- new_prmt[[2]]
        prmt$mu<- new_prmt[[3]]
        prmt$sigma <- new_prmt[[4]]
        
        # New log likelihood
        mixture_model <- prmt$lambda * dexp(x, rate = prmt$theta, log = FALSE) + 
          (1 - prmt$lambda) * dnorm(x, mean = prmt$mu, sd = prmt$sigma)
        loglik = sum(log10(mixture_model))
        
        eps <- (loglik - loglik_old)^2
        if (is.na(eps)) {
          eps <- 10
        } else {
          loglik_old <- loglik
        }
        iter <- iter + 1
      }
      
      # Calculate d_i for each cell and update D table
      exp_distr = prmt$lambda * dexp(x, rate = prmt$theta, log = FALSE)
      norm_distr = (1 - prmt$lambda) * dnorm(x, mean = prmt$mu, sd = prmt$sigma)
      D[gene, ] = exp_distr / (exp_distr + norm_distr)
      
      parameters[gene, "lambda"] <- prmt$lambda
      parameters[gene, "theta"] <- prmt$theta
      parameters[gene, "mu"] <- prmt$mu
      parameters[gene, "sigma"] <- prmt$sigma
    }
  }
  result <- list("D"=D, "parameters"=parameters)
  return(result)
}

dropout_correction <- function(gem) {
  # Calculate parameters and dropout probability with exponential distr
  exp_mix_model <- EM_for_dropout_estimation(gem)
  
  D <- 1 - exp_mix_model$D
  
  # normalize gem
  gem_norm <- t(apply(gem, 1, normalized_exp_pctl))
  
  # Get gene presence matrix corrected for dropouts
  gpm <- gem
  for (gene in rownames(gem_norm)) {
    # Also added na.rm=TRUE
    drop_cand <- which(gem_norm[gene, ] <= quantile(gem_norm[gene, ], 0.25, na.rm = TRUE))
    for (c in drop_cand) {
      gpm[gene, c] <- D[gene, c]
    }
  }
  
  return(gpm)
}
