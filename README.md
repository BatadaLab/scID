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





