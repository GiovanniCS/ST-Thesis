# Application of too-many-cells clustering algorithm
In this experiment the divisive hierarchical clustering algorithm developed for too-many-cells is applied to a spatial transciptomic dataset.

The dataset is the slice "anterior1" provided in the stxBrain data package of Seurat.

Workflow:

R
````R
library(Seurat)
library(SeuratData)
library(TooManyCellsR)

brain <- LoadData("stxBrain", type = "anterior1")
counts <- brain@assays$Spatial@counts
writeMatrixFiles(counts) #TooManyCellsR packages
#locate genes.tsv, barcodes.tsv and matrix.mtx
````
bash
```bash
cd /Users/giovanni/Desktop/ST-Thesis/t-m-cOnMouseBrain/input

sudo docker run -it --rm -v "/Users/giovanni/Desktop/ST-Thesis/t-m-cOnMouseBrain/input:/home" gregoryschwartz/too-many-cells:0.2.2.0 make-tree --matrix-path /home --output /home/out > clusters.csv
```

Read clusters.csv into R
```R
clusters <- read.csv("clusters.csv",header = FALSE)
new_clusters <- factor(clusters$V2,labels = clusters$V1)
names(new_clusters) <- clusters$V1
brain@active.ident <- new_clusters
SpatialDimPlot(brain, label = TRUE, label.size = 3) + NoLegend()
```