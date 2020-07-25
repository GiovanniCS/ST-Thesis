# Seurate distance measure extension
In this experiment the distance measure computation of Seurat is bypassed. Intead, a distance matrix is provided to the FindNeighbors function. This matrix is the same that FindNeighbors would compute with the addition of fisical distance between spots on the chip.

Workflow:
````R
library(Seurat)
library(RColorBrewer)

mouse <- Load10X_Spatial(data.dir = "./filtered_feature_bc_matrix")
mouse <- SCTransform(mouse, assay = "Spatial", verbose = FALSE)
mouse <- RunPCA(mouse, assay = "SCT", verbose = FALSE)
m = mouse@reductions[["pca"]]@cell.embeddings
distPCA = dist(m,diag=TRUE)  
distPCA = as.matrix(distPCA)
coord <- GetTissueCoordinates(mouse)
distcoord = dist(coord,diag=TRUE)
distcoord = as.matrix(distcoord)
distcoord = log(distcoord+1)
finaldist = distPCA + distcoord
finaldist = as.dist(finaldist)
neighbors = FindNeighbors(finaldist)
neighbors = list(neighbors_nn=neighbors$nn,neighbors_snn=neighbors$snn)
mouse@graphs = neighbors
mouse <- FindClusters(mouse, verbose = FALSE,graph.name = "neighbors_snn")

palette = colorRampPalette(brewer.pal(11,name="BrBG"))
clusters = length(levels(mouse@active.ident))
SpatialDimPlot(mouse, label = TRUE, label.size = 3,cols=palette(clusters))
````

Default Seurat clustering
![](https://user-images.githubusercontent.com/25981629/88466213-4e43c000-ceca-11ea-9ec8-a444ab91e48b.png)

New Seurat clustering with fisical distance added
![](https://user-images.githubusercontent.com/25981629/88466220-5a2f8200-ceca-11ea-96e3-402bdedbfd2f.png)



## Other datasets
### Human heart
Default Seurat clustering
![](https://user-images.githubusercontent.com/25981629/88466563-f6a75380-cecd-11ea-8814-04e267ea8d3b.png)
New Seurat clustering with fisical distance added
![](https://user-images.githubusercontent.com/25981629/88466566-06269c80-cece-11ea-90d5-999fa5863be0.png)

### Human limph node
Default Seurat clustering
![](https://user-images.githubusercontent.com/25981629/88466602-530a7300-cece-11ea-9f63-202a2f94838f.png)
New Seurat clustering with fisical distance added
![](https://user-images.githubusercontent.com/25981629/88466624-94028780-cece-11ea-854d-794c7500a8ea.png)

### Mouse kidney
Default Seurat clustering
![](https://user-images.githubusercontent.com/25981629/88466699-4d615d00-cecf-11ea-8722-64c73cd5b046.png)
New Seurat clustering with fisical distance added
![](https://user-images.githubusercontent.com/25981629/88466696-44708b80-cecf-11ea-82e9-0cd6749ac378.png)