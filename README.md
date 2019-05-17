# scID

The power of single cell RNA sequencing (scRNA-seq) stems from its ability to uncover cell type-dependent phenotypes, which rests on the accuracy of cell type identification. However, resolving cell types within and, thus, comparison of scRNA-seq data across conditions is challenging due to technical factors such as sparsity, low number of cells and batch effect. To address these challenges we developed scID (Single Cell IDentification), which uses the framework of Fisher's Linear Discriminant Analysis to identify transcriptionally related cell types between scRNA-seq datasets. We demonstrate in the preprint ["Mapping transcriptionally equivalent populations across single cell RNA-seq datasets"](https://www.biorxiv.org/content/10.1101/470203v1) the accuracy and performance of scID relative to existing methods on several published datasets. By increasing power to identify transcriptionally similar cell types across datasets, scID enhances investigator's ability to extract biological insights from scRNA-seq data.

## Installation
scID can be installed using the devtools R package:
```
devtools::install_github("BatadaLab/scID")
```




