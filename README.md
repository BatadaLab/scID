# scID

The power of single cell RNA sequencing (scRNA-seq) stems from its ability to uncover cell type-dependent phenotypes, which rests on the accuracy of cell type identification. However, resolving cell types within and, thus, comparison of scRNA-seq data across conditions is challenging due to technical factors such as sparsity, low number of cells and batch effect. To address these challenges we developed scID (Single Cell IDentification), which uses the framework of Fisher's Linear Discriminant Analysis to identify transcriptionally related cell types between scRNA-seq datasets. We demonstrate the accuracy and performance of scID relative to existing methods on several published datasets. By increasing power to identify transcriptionally similar cell types across datasets, scID enhances investigator's ability to extract biological insights from scRNA-seq data.

## How scID works
scID uses Linear Discriminant Analysis approach described in the preprint "" link

## Installation
scID can be installed using the devtools R package:
```
devtools::install_github("BatadaLab/scID")
```

## Usage
```
scid_match_cells(target_gem=T, reference_gem=R, reference_clusters=L, ...)
```

### Input
scID requires the following inputs:
1. A data frame of TPM-normalised gene (rows) expression in cells (columns) of the target data (R)
2. A data frame of TPM-normalised gene (rows) expression in cells (columns) of the reference data (T) and a list of cluster IDs for the reference cells or alternatively (L), 
3. A data frame of signature genes per reference cluster. This has to have columns named "gene" and "cluster" for the gene symbols and cluster IDs for which each gene is a marker.

### Output
1. A list of labels for the target cells
2. A data frame of signature genes per reference cluster.

## Tutorial
This tutorial is an example pipeline for mapping across two 10X datasets of E18 mouse brain cells and nuclei from cortex, hippocampus and subverticular zone. The raw counts for the reference and the target data can be found [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/neuron_9k) and [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/nuclei_900) respectively. However, here we provide the TPM normalized values to speed-up preprocessing.

The reference cells can be grouped into 15 clusters as shown in the tSNE plot
![](https://github.com/BatadaLab/scID/blob/master/ExampleData/figures/Reference_tSNE.png)


After loading the libraries we read the files
```
library(scID)

target_gem <- readRDS(file="~/scID/ExampleData/target_gem.rds")

reference_gem <- readRDS(file="~/scID/ExampleData/reference_gem.rds")
reference_clusters <- readRDS(file="~/scID/ExampleData/reference_clusters.rds")

scID_res <- scid_match_cells(target_gem = target_gem, reference_gem = reference_gem, 
                             reference_clusters = reference_clusters, logFC = 0.5, likelihood_threshold = 0.95)
```

Alternatively, scID can take a data frame of signature genes per cluster without reference cells. This could also be curated lists of markers. 
```
markers <- readRDS(file="~/scID/ExampleData/markers.rds")

scID_res <- scid_match_cells(target_gem = target_gem, markers = markers, logFC = 0.5, likelihood_threshold = 0.95)
```

The next heatmap shows the average expression of each markers' list in each of the reference clusters. Each row represents a markers' list and each column a cluster of cells.
```
make_heatmap(gem = reference_gem, labels = reference_labels, markers = markers)
```
![](https://github.com/BatadaLab/scID/blob/master/ExampleData/figures/Reference_heatmap.png)

The respective heatmap of target nuclei data grouped by scID can is shown below
```
make_heatmap(gem = target_gem, labels = scID_res$labels, markers = markers)
```
![](https://github.com/BatadaLab/scID/blob/master/ExampleData/figures/Target_heatmap.png)

Here we can see that the gene expression pattern in the target data is very similar to the one in the reference data. The gray column shows that no target cell were assigned to the reference clusters 5 and 13.




