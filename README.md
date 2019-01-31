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

scID requires the following inputs:
1. A data frame of TPM-normalised gene (rows) expression in cells (columns) of the target data
2. A data frame of TPM-normalised gene (rows) expression in cells (columns) of the reference data and a list of cluster IDs for the reference cells or alternatively, 
3. A data frame of markers per reference cluster. This has to have columns named "gene" and "cluster" for the gene symbols and cluster IDs for which each gene is a marker.

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

Here we can see that the gene expression pattern in the target data is very similar to the one in the reference data. The gray column shows that no target cell was assigned to the reference "Cluster 3".




