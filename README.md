# scID

The power of single cell RNA sequencing (scRNA-seq) stems from its ability to uncover cell type-dependent phenotypes, which rests on the accuracy of cell type identification. However, resolving cell types within and, thus, comparison of scRNA-seq data across conditions is challenging due to technical factors such as sparsity, low number of cells and batch effect. To address these challenges we developed scID (Single Cell IDentification), which uses the framework of Fisher's Linear Discriminant Analysis to identify transcriptionally related cell types between scRNA-seq datasets. We demonstrate in the preprint ["Mapping transcriptionally equivalent populations across single cell RNA-seq datasets"](https://www.biorxiv.org/content/10.1101/470203v1) the accuracy and performance of scID relative to existing methods on several published datasets. By increasing power to identify transcriptionally similar cell types across datasets, scID enhances investigator's ability to extract biological insights from scRNA-seq data.

scID classifies cells of a given target dataset based on their transcriptional similarity to given reference clusters in 4 steps. As a first step, scID extracts cluster-specific gene sets from the reference data and calculates weights (based on Fisher's Linear Discriminant Analysis) that represent their discriminative power to identify the cluster of interest. Next, scID scores all target cells based on the expression of the cluster-specific gene sets and, finally, identifies equivalent target cells by fitting a mixture of Gaussian distributions. 

![](https://github.com/BatadaLab/scID/blob/master/scID_pipeline.png)


## Installation
scID can be installed using the devtools R package:
```
devtools::install_github("BatadaLab/scID")
```




