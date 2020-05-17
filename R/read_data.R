#' Function to read Gene Expression Matrix from file
#' 
#' @param filename directory and name of input file
#' @param header specify if first line of file is the field names
#' 
#' @return Data frame of GEM with genes in rows and cells in columns
#' @export
loadfast <- function(filename, header = T) {
  
  df <- data.table::fread(filename, header=header, showProgress=TRUE, data.table=FALSE)
  rownames(df) <- toupper(df[,1])
  df <- df[,-1]
  
  df
}

#' Function to read markers from .gmt file format (Genomic Cytometry)
#' 
#' @param filename .gmt file with markers
#' 
#' @return Data frame of cluster specific genes per celltype
#' @export
gmt_to_markers <- function(gmt_file) {
  
  data <- qusage::read.gmt(gmt_file)
  celltypes <- names(data)
  markers <- data.frame(gene = NA, cluster = NA, avg_logFC = NA)
  for (ct in celltypes) {
    
    m <- data.frame(gene = data[[ct]], cluster = ct, avg_logFC = 1)
    markers <- rbind(markers, m)
  }
  markers <- markers[complete.cases(markers), ]
  
  markers
}