# scID

The power of single cell RNA sequencing (scRNA-seq) stems from its ability to uncover cell type-dependent phenotypes, which rests on the accuracy of cell type identification. However, resolving cell types within and, thus, comparison of scRNA-seq data across conditions is challenging due to technical factors such as sparsity, low number of cells and batch effect. To address these challenges we developed scID (Single Cell IDentification), which uses the framework of Fisher's Linear Discriminant Analysis to identify transcriptionally related cell types between scRNA-seq datasets. We demonstrate in the preprint ["Mapping transcriptionally equivalent populations across single cell RNA-seq datasets"](https://www.biorxiv.org/content/10.1101/470203v1) the accuracy and performance of scID relative to existing methods on several published datasets. By increasing power to identify transcriptionally similar cell types across datasets, scID enhances investigator's ability to extract biological insights from scRNA-seq data.

## Installation
scID can be installed using the devtools R package:
```
devtools::install_github("BatadaLab/scID")
```

## Usage

There are three ways to use scID. 

### Usage 1: Canonical usage
Given two datasets of single-cell RNA-seq gene expression for which cell grouping for one the datasets (reference) is known, scID seeks to find transcriptionally equivalent groups of cells for the second dataset (target).
```
scID_output <- scID::scid_multiclass(target_gem, reference_gem, reference_clusters, ...)
```

#### Input
1. ```target_gem``` An nxm data frame of n genes (rows) in m cells (columns) of the dataset with unknown grouping, where each entry is library-depth or column normalized gene expression. Cell names are expected to be unique.
2. ```reference_gem``` An NxM data frame of N genes (rows) in M cells (columns) of the dataset with known grouping, where each entry is library-depth or column normalized gene expression. 
3. ```reference_clusters``` A list of cluster labels for the reference cells.

#### Output

scID_output is a list of two objects 

1. ```scID_output$labels``` A named list of cluster labels for the target cells

2. ```scID_output$markers``` A data frame of signature genes extracted from the reference clusters. 

### Usage 2: Single reference but multiple targets (T1, T2)

* Step 1: Extract markers from reference clusters
```
markers_generated_by_scID <- scID::find_markers(reference_gem, reference_clusters, logFC)
```
This step can be skipped when the user has own method for extracting markers. 

* Step 2: Find transctiptionally equivalent cells in target datasets
```
scID_output_T1 <- scID:scid_multiclass(T1, markers_generated_by_scID, ...)

scID_output_T2 <- scID::scid_multiclass(T2, markers_generated_by_scID, ...)
```
### Usage 3: User-specified cluster gene signatures
A pre-computed set of markers can be given as input by the user alternatively. The markers object has to be a data frame with genes and cluster ID in columns as in [this](https://github.com/BatadaLab/scID/blob/master/ExampleData/markers.rds) example file.
```
scID_output <- scID::scid_multiclass(T, markers_generated_by_user, ...)
```


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

scID_output <- scid_multiclass(target_gem = target_gem, reference_gem = reference_gem, 
                               reference_clusters = reference_clusters, logFC = 0.5, likelihood_threshold = 0.95)
```

Alternatively, scID can take a data frame of signature genes per cluster without reference cells. This could also be curated lists of markers. 
```
markers <- readRDS(file="~/scID/ExampleData/markers.rds")

scID_output <- scid_multiclass(target_gem = target_gem, markers = markers, logFC = 0.5, likelihood_threshold = 0.95)
```

The next heatmap shows the average expression of each markers' list in each of the reference clusters. Each row represents a markers' list and each column a cluster of cells.
```
make_heatmap(gem = reference_gem, labels = reference_clusters, markers = markers)
```
![](https://github.com/BatadaLab/scID/blob/master/ExampleData/figures/Reference_heatmap.png)

The respective heatmap of target nuclei data grouped by scID can is shown below
```
make_heatmap(gem = target_gem, labels = scID_output$labels, markers = markers)
```
![](https://github.com/BatadaLab/scID/blob/master/ExampleData/figures/Target_heatmap.png)

Here we can see that the gene expression pattern in the target data is very similar to the one in the reference data. The gray column shows that no target cell were assigned to the reference clusters 5 and 13.




