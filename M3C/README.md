# M3C package application to mouse coronal brain section
In this experiment the M3C R package is used to cluster the dataset with PAM algorithm.  
M3C estimate the most suitable number of clusters to be discovered expliting the concept of consensus clustering based on subsets of the orignal observation set.  
The clustering algo. available are PAM (default), K-means and hierarchical spectral clustering.  

Workflow:
````R
library(Seurat)
library(M3C)

mouse <- Load10X_Spatial(data.dir = "./Datasets/MouseBrainCoronal/filtered_feature_bc_matrix")
mouse <- SCTransform(mouse, assay = "Spatial", verbose = FALSE)
mydatamouse = matrix(mouse@assays[["SCT"]]@counts,18768,2702)
res <- M3C(mydatamouse, removeplots = TRUE, iters=25, objective='PAC', fsize=8, lthick=1, dotsize=1.25)
````

Expert annotation of the coronal mouse brain section
![](https://user-images.githubusercontent.com/25981629/88228893-248f5c80-cc70-11ea-91de-61aa89601e66.png)
