#' Function to read Gene Expression Matrix from file
#' @param filename directory and name of input file
#' @param header specify if first line of file is the field names
#' @return Data frame of GEM with genes in rows and cells in columns
#' @export
loadfast <- function(filename, header=T) {
  ## why not use data.table? must faster for gem
  require(data.table) # https://www.analyticsvidhya.com/blog/2016/05/data-table-data-frame-work-large-data-sets/
  df=fread(filename, header=header, showProgress=TRUE, data.table=FALSE)
  rownames(df)=toupper(df[,1])
  df=df[,-1]
  return(df)
}