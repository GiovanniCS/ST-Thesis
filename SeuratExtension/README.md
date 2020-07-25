# Seurate distance measure extension
In this experiment the distance measure computation of Seurat is bypassed. Intead, a distance matrix is provided to the FindNeighbors function. This matrix is the same that FindNeighbors would compute with the addition of fisical distance between spots on the chip.

Workflow:
````R
library(Seurat)
library(RColorBrewer)

mouse <- Load10X_Spatial(data.dir = "./Datasets/MouseBrainCoronal/filtered_feature_bc_matrix")
mouse <- SCTransform(mouse, assay = "Spatial", verbose = FALSE)
mouse <- RunPCA(mouse, assay = "SCT", verbose = FALSE)
m = mouse@reductions[["pca"]]@cell.embeddings
distPCA = dist(m,diag=TRUE)  
distPCA = as.matrix(distPCA)
coord <- GetTissueCoordinates(mouse)
distcoord = dist(coord,diag=TRUE)
distcoord = as.matrix(distcoord)
finaldist = distPCA + distcoord
mouse <- FindNeighbors(mouse, distance.matric = finaldist)
mouse <- FindClusters(mouse, verbose = FALSE)

palette = colorRampPalette(brewer.pal(11,name="BrBG"))
clusters = length(levels(mouse@active.ident))
SpatialDimPlot(mouse, label = TRUE, label.size = 3,cols=palette(clusters))
````

Default Seurat clustering
![](https://user-images.githubusercontent.com/25981629/88464572-349b7c00-cebc-11ea-9845-054ea8ef343a.png)

New Seurat clustering with fisical distance added
![](https://user-images.githubusercontent.com/25981629/88464570-2baaaa80-cebc-11ea-9e8f-f8779a03454c.png)