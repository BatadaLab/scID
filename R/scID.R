#' Main function to get Gene expression matrix and signature genes and return matches and scores
#' @param signature_file (optional) filename of the signature genes
#' @param signature_genes (optional) a list of the signature genes
#' either signature_file or signature_genes must be given as input
#' @param gem_file (optional) filename of the gene expression matrix
#' @param scData (optional) Data frame of gene expression (rows) per cell (columns)
#' either gem_file or scData must be given as input
#' @param positive_markers list of genes that should be present in confident IN cells (optional)
#' @param negative_markers list of genes that should not be present in any IN cell (optional)
#' @param contamination Percentage of accepted cells that belong to the common area between 
#' IN and OUT population distributions
#' @param species Choose between "human" and "mouse" to use a predefined list of housekeeping
#' genes for cell quality correction
#' @param hk_genes List of houskeeping genes in case species not defined
#' @param sort.signature Logical to choose to skip gene sorting and give equal weights to all genes
#' (default is TRUE for sorting cells)
#' @param do.imputation Logical to choose not to correct for dropouts (default is TRUE for doing imputation)
#' @return list of names of matching cells (IN-population) and matching score of every cell of the dataset
#' @return list of weight per gene
#' @export
scid_match_cells <- function(signature_file=NULL, gem_file=NULL, scData=NULL, signature_genes=NULL, 
                             positive_markers=NULL, negative_markers=NULL, contamination=0,  
                             species = "human", hk_genes = NULL, sort.signature = TRUE, do.imputation = TRUE) {
  
  source("dropout_correction.R")
  source("read_data.R")
  source("rank_signature_genes.R")
  source("calculate_score.R")
  
  ## embed the read.delim.with.error.handling() inside this function to make it private
  
  # Read housekeeping genes
  if (!is.null(hk_genes)) {
    hk_genes <- toupper(hk_genes)
  } else if (species=="human") {
    hk_genes <- toupper(read.table("../data/housekeeping_hsa.txt", stringsAsFactors = F)$V1)
  } else if(species=="mouse") {
    hk_genes <- toupper(read.table("../data/housekeeping_mmu.txt", stringsAsFactors = F)$V1)
  } else {
    print("Species not known, please choose between 'human' and 'mouse' or provide a list of housekeeping genes")
  }
  
  #----------------------------------------------------------------------------------------------------
  # Read data
  if (is.null(scData)) {
    scData <- loadfast(gem_file)
  }
  if (is.null(signature_genes)) {
    signature_genes <- toupper(read.table(signature_file, stringsAsFactors = FALSE)$V1)
  }
  
  # Make row and column names capital
  rownames(scData) <- toupper(rownames(scData))
  #colnames(scData) <- toupper(colnames(scData))
  
  # Keep only genes in scData
  signature_genes <- intersect(signature_genes, rownames(scData))
  
  #----------------------------------------------------------------------------------------------------
  # Find score for signature genes
  # Calculate specificity
  if (sort.signature) {
    weights <- weight_signature(gem=scData, signature = signature_genes, positive_markers, negative_markers)
  } else {
    weights <- rep(1, length(signature_genes))
    names(weights) <- signature_genes
  }
  if (typeof(weights)=="character") {
    print(weights)
  } else {
    # Log transform data (log(gem+1.01))
    gem <- scData[names(which(weights > 0)), ]
    gem <- log2(apply(gem, c(1,2), function(x) 2^x-1) + 1.01)
    
    # Calculate dropout probability
    if (do.imputation) {
      gpm <- dropout_correction(gem)
    } else {
      gpm <- t(apply(gem, 1, normalized_exp_pctl))
    }
    # Remove NA values
    na.values <- which(apply(gpm, 1, function(x){any(is.na(x))}))
    if (length(na.values > 0)) {
      gpm <- gpm[-na.values, ]
    }
    # Calculate matching score
    genes <- intersect(names(weights), rownames(gpm))
    weights <- weights[genes]
    gpm <- gpm[genes, ]
    
    matching_score <- colSums(na.omit(weights*gpm))/sum(na.omit(weights))
    
    adjusted_score <- adjust_score(scData = scData, matching_score = matching_score, hk_genes = hk_genes)
    
    populations <- final_populations(adjusted_score, contamination = contamination)
    
    if (!sort.signature | (length(populations$IN) == 0)) {
      print(paste("Found", length(populations$IN), "cells matching"))
      return(list(matches=populations$IN, matchingScore=adjusted_score, geneWeights=weights))
    } else {
      confident_IN <- populations$IN
      if (length(populations$OUT) > 0) {
        confident_OUT <- populations$OUT
      } else {
        rest_cells <- setdiff(colnames(scData), confident_IN)
        confident_OUT <- rest_cells[sample(length(rest_cells), 0.2*length(rest_cells))]
      }
      
      # Re-weight signature with better IN and OUT cells
      snr <- calculate_snr(na.omit(scData[genes, ]), true_cells = confident_IN, false_cells = confident_OUT)
      weights <- snr/max(snr)
      
      matching_score <- colSums(na.omit(weights*gpm))/sum(na.omit(weights))
      adjusted_score <- adjust_score(scData = scData, matching_score = matching_score, hk_genes = hk_genes)
      
      populations <- final_populations(adjusted_score, contamination = contamination)
      
      print(paste("Found", length(populations$IN), "cells matching"))
      
      return(list(matches=populations$IN, matchingScore=adjusted_score, geneWeights=weights))
    }
  }
}

