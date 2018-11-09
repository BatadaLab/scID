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

scID input should be a TPM-normalised Gene Expression Matrix and a list of genes that describe the population of interest. Both datasets can be read from files or given as ready R objects (matrix/data.frame/list). 

## Tutorial
This tutorial is an example pipeline for mapping across two 10X datasets of E18 mouse brain cells and nuclei from cortex, hippocampus and subverticular zone. The raw counts for the reference and the target data can be found here[Reference](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/neuron_9k) and here[target](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/nuclei_900) respectively. 


* For providing datasets already loaded on R environment:
```
scID_output <- scid_match_cells(signature_genes = signature_list, gem = GEM_dataframe)
```

* Output consists of three objects; the selected matching cells, a list of matching scores for all cells of the dataset, a list of gene weights generated during the sorting step.
```
matching_cells <- scID_output$matches
matching_scores <- scID_output$matchingScore
gene_weights <- scID_output$gene.weights
```







