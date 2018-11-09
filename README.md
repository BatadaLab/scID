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
library(pheatmap)


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
![](https://github.com/BatadaLab/scID/blob/master/ExampleData/figures/Reference_tSNE.png)

Next, we find positive markers of these 15 clusters using MAST.
```
markers <- FindAllMarkers(sobj_ref, only.pos = TRUE, test.use = "MAST", logfc.threshold = 0.5)
```
The next heatmap shows the average expression of each markers' list in each of the reference clusters. Each row represents a markers' list and each column a cluster of cells.
```
gem_avg_ref <- data.frame(matrix(NA, length(unique(sobj_ref@ident)), length(unique(sobj_ref@ident))), 
                          row.names = paste("Cluster", unique(sobj_ref@ident), "geneset", sep = "_"))
colnames(gem_avg_ref) <- paste("Cluster", unique(sobj_ref@ident), sep = "_")
for (i in 1:length(unique(sobj_ref@ident))) {
  cells <- WhichCells(sobj_ref, i-1)
  if (length(unique(sobj_ref@ident)) > 1) {
    avg_exp <- rowMeans(reference_gem[markers$gene, cells])
  } else if (length(cells) == 1) {
    avg_exp <- gem[markers$gene, cells]
    names(avg_exp) <- markers$gene
  } else {
    next
  }
  for (j in 1:length(unique(sobj_ref@ident))) {
    gem_avg_ref[j,i] <- mean(na.omit(avg_exp[markers$gene[which(markers$cluster == j-1)]]))
  }
}

pheatmap(gem_avg_ref, border="white", color = colorspace::diverge_hsv(50), cluster_rows = F,
         cluster_cols = F, show_colnames = T, fontsize_row = 10, border_color = F,show_rownames = T,
         scale = "row")
```
![](https://github.com/BatadaLab/scID/blob/master/ExampleData/figures/Reference_heatmap.png)

Now that we have extracted features of the reference clusters we can map the target nuclei data to these clusters with scID. 
For multiclass mapping we propose the following approach:

* **Step 1:** We run scID for each of the sets of features and store the scores of the matching cells in a table.
```
scores <- data.frame(matrix(NA, length(unique(sobj_9k@ident)), ncol(target_gem)), row.names = unique(sobj_9k@ident))
colnames(scores) <- colnames(target_gem)
for (celltype in unique(sobj_9k@ident)) {
  signature <- markers$gene[which(markers$cluster == celltype)]
  res <- scid_match_cells(gem = target_gem, signature_genes = signature, contamination = 0.05)
  scores[celltype, res$matches] <- scale(res$matchingScore[res$matches])
}
```

* **Step 2:** Now we need to resolve any conflicts. In datasets that are sparse and/or have very similar subtypes of cells, it is common for a cell to be selected as matching to more that one signatures. To resolve this issue, we need to make the cell’s scores for the different sets of features for which it was a match comparable. Thus why in the previous step, the scores of matching cells for each set of features are centered and we resolve conflicting classes by assigning each cell to the class with the highest score. Cells that were not selected by any of the gene lists are marked as "unassigned".

```
scID_labels <- c()
for (cell in colnames(scores)) {
  if (all(is.na(scores[, cell]))) {
    matching_type <- "unassigned"
  } else {
    matching_type <- rownames(scores)[which(scores[, cell] == max(scores[, cell], na.rm = T))]
  }
  scID_labels <- c(scID_labels, matching_type)
}
names(scID_labels) <- colnames(scores)
```

Now we can visualize the average expression of these gene lists in each group of target cells as clustered by scID. 

```
markers_filt <- markers[which(markers$gene %in% rownames(target_gem)), ]
gem_avg_targ <- data.frame(matrix(NA, length(unique(sobj_ref@ident)), length(unique(sobj_ref@ident))),
                      row.names = paste("Cluster", unique(sobj_ref@ident), "geneset", sep = "_"))
colnames(gem_avg_targ) <- paste("Cluster", unique(sobj_ref@ident), sep = "_")

for (i in 1:length(unique(sobj_ref@ident))) {
  cells <- names(scID_labels)[which(scID_labels == i-1)]
  if (length(cells) > 1) {
    avg_exp <- rowMeans(target_gem[markers_filt$gene, cells])
  } else if (length(cells) == 1) {
    avg_exp <- target_gem[markers_filt$gene, cells]
    names(avg_exp) <- markers_filt$gene
  } else {
    next
  }
  for (j in 1:length(unique(sobj_ref@ident))) {
    gem_avg_targ[j,i] <- mean(na.omit(avg_exp[markers_filt$gene[which(markers_filt$cluster == j-1)]]))
  }
}

pheatmap(gem_avg_targ, border="white",color = colorspace::diverge_hsv(50), cluster_rows = F,
         cluster_cols = F, show_colnames = T, fontsize_row = 10, border_color = F,show_rownames = T,
         scale = "row")
```
![](https://github.com/BatadaLab/scID/blob/master/ExampleData/figures/Target_heatmap.png)

Here we can see that the gene expression pattern in the target data is very similar to the one in the reference data. The gray column shows that no target cell was assigned to the reference "Cluster 3".




