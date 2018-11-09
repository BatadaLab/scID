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
This tutorial is an example pipeline for mapping across two 10X datasets of E18 mouse brain cells and nuclei from cortex, hippocampus and subverticular zone. The raw counts for the reference and the target data can be found [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/neuron_9k) and [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/nuclei_900) respectively. However, here we provide the TPM normalized values to speed-up preprocessing.

After loading the libraries we read the files
```
library(scID)
library(Seurat)

reference_gem <- readRDS(file="~/scID/ExampleData/reference_gem.rds")
target_gem <- readRDS(file="~/scID/ExampleData/target_gem.rds")
```

scID can take a list of features without reference cells, e.g. curated lists of markers. However, often such information is not available, thus, we can extract the features from a clustered reference dataset.

Here, we will cluster the reference data using [Seurat](https://satijalab.org/seurat/) and extract markers for each cluster.
```
sobj_ref <- CreateSeuratObject(raw.data = reference_gem)
sobj_ref <- NormalizeData(sobj_ref)
sobj_ref <- ScaleData(sobj_ref)
sobj_ref <- FindVariableGenes(sobj_ref)
sobj_ref <- RunPCA(sobj_ref,  do.print = FALSE)
sobj_ref <- ProjectPCA(sobj_ref)
sobj_ref <- FindClusters(sobj_ref, dims.use = 1:5, save.SNN = F)
sobj_ref <- RunTSNE(sobj_ref, dims.use = 1:5, do.fast = T)

TSNEPlot(sobj_ref, do.label = T, pt.size = 0.1, no.axes = T, no.legend = T)
```
This results in 15 clusters as shown in the tSNE plot
![tSNE]("ExampleData/figures/Reference_tSNE.png")








