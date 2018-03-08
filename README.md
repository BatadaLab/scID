# scID

scID is a method for geneset-guided identification of cell types at the level of individual cells, using single cell RNA-seq data.
Given a gene expression matrix and a list of marker genes that define the population of interest, scID returns a matching score for each cell and the names of the identified matching cells.

One of the advantages of scID is gene prioritisation; the signature genes are sorted by importance in distinguishing the cell type of interest from the background cells. 
This is tailored to every specific experiment, accounting for sequencing depth and dropouts, as well as diversity of populations in the mixture. 
Additionally, it can handle contaminated signatures with false markers, so there is no requirement for a very confident and "clean" signature.
Finally scID accounts and corrects for cell quality differences to make score comparable between cells of different quality.

## Installation
scID can be installed using the devtools R package:

```
devtools::install_github("BatadaLab/scID")
```

## Input parameters

scID input should be a TPM-normalised Gene Expression Matrix and a list of genes that describe the population of interest. Both datasets can be read from files or given as ready R objects (matrix/data.frame/list). Additionally, the user can specify a subset of the gene signature that required to be expressed in all cells of the IN-population (positive markers) and a set of negative markers that should not be expressed in any cells of the IN-population. Both sets of genes can help scID more accuratelly prioritize the signature genes when the background is very similar to the target population, but are not necessary. 
To account for cell quality differences, scID uses the expression of housekeeping genes in the cells. We currently provide lists of housekeeping genes for human and mouse datasets by specifying species as input (default is human). When gene expression comes from a different species you need to also provide the list of housekeeping genes for your orgnaism.

## Example code





