# Space-aware clustering 
In this out we try to apply the findings described in [Efficient spatial segmentation ...](https://pdfs.semanticscholar.org/51f9/b2095466e70edf3add5993c0c1c7800e79a5.pdf?_ga=2.236486721.1509296317.1594040010-1511101378.1594040010) to the dataset 1. [Here](https://www.youtube.com/watch?v=_3U2Elt5CTI) is the presentation of this paper by the author.  
Dataset 1: slice "anterior1" provided in the stxBrain data package of Seurat.

Although this clustering algorithm was developed for Mass Spectrometry Imaging data, it could theoretically be applied to other kind of spatially resolved data ( such as spatial transcriptomics) as suggested in the paper --> "Our spatial segmentation methods can be applied for segmenting other hyper-spectral or multi-channel data"

One open source implementation of this method is the R package [Cardinal](https://www.bioconductor.org/packages/release/bioc/html/Cardinal.html), focused on the analysis of mass spectrometry imaging datasets.

Workflow:

````R
library(Seurat)
library(SeuratData)
library(Cardinal)

brain <- LoadData("stxBrain", type = "anterior1")
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

#the code below simulate a mass spectrometry imaging dataset with spatial transcriptomics data of brain dataset
m = matrix(brain@assays$SCT@counts,17668,2696)
coord <- GetTissueCoordinates(brain)
run <- factor(rep("run0", nrow(coord)))
s <- simulateSpectrum(n=1, peaks=10, from=100, to=150000)
mz = s$mz[1:17668]
fdata <- MassDataFrame(mz)
pdata <- PositionDataFrame(run=run, coord=coord)
out <- MSImagingExperiment(imageData=m,featureData=fdata,pixelData=pdata)

#Space-aware clustering based on k-means clustering --> original solution of Alexandrov
skm = spatialKMeans(out,r=5,k=15,method="adaptive")
clusters = skm@resultData@listData[[1]][["cluster"]]
names(clusters) = names(brain@active.ident)
brain@active.ident <- clusters
SpatialDimPlot(brain, label = TRUE, label.size = 3) + NoLegend()
````

## Seurat reference:
![](https://user-images.githubusercontent.com/25981629/86606773-b9276880-bfa8-11ea-85e8-db9c498b85eb.png)

#Difference between gaussian and adaptive weights:
![](https://user-images.githubusercontent.com/25981629/86608264-926a3180-bfaa-11ea-87c0-f7f95996f2d4.png)

## Results with adaptive weights
![](https://user-images.githubusercontent.com/25981629/86606428-43230180-bfa8-11ea-9764-8b86b38489f3.png)
![](https://user-images.githubusercontent.com/25981629/86604816-271e6080-bfa6-11ea-9fb2-93f2167848c8.png)
![](https://user-images.githubusercontent.com/25981629/86603581-7368a100-bfa4-11ea-9dcf-3f8a1bc94773.png)
![](https://user-images.githubusercontent.com/25981629/86601418-92196880-bfa1-11ea-8370-7564a636ce1f.png)

## Results with gaussian weights
![](https://user-images.githubusercontent.com/25981629/86607746-f3453a00-bfa9-11ea-95ae-c1d5cb97c447.png)
![](https://user-images.githubusercontent.com/25981629/86608839-597e8c80-bfab-11ea-86de-56202dab1761.png)
![](https://user-images.githubusercontent.com/25981629/86614523-7028e180-bfb3-11ea-981d-df0a6964d92f.png)
![](https://user-images.githubusercontent.com/25981629/86616136-b3844f80-bfb5-11ea-8a64-dda9ef8cae93.png)

````R
#Alternative space-aware clustering proposed by Cardinal
ssc <- spatialShrunkenCentroids(out, r=2, k=14, s=c(0,3,6,9))
clusters = ssc@resultData@listData[[1]][["cluster"]]
names(clusters) = names(brain@active.ident)
brain@active.ident <- clusters
SpatialDimPlot(brain, label = TRUE, label.size = 3) + NoLegend()
````
## Results with gaussian weigths, s = 0,3,6 and r = 1,2 (s = sparsity thresholding parameter by which to shrink the t-statistics)
![](https://user-images.githubusercontent.com/25981629/86618282-eda32080-bfb8-11ea-8ffe-79e5015c0cf0.png)
![](https://user-images.githubusercontent.com/25981629/86618375-19260b00-bfb9-11ea-9395-4167745258f0.png)
![](https://user-images.githubusercontent.com/25981629/86618484-4a9ed680-bfb9-11ea-9469-29a24eef1404.png)
![](https://user-images.githubusercontent.com/25981629/86619916-b3874e00-bfbb-11ea-996b-34c0b6cd18cb.png)
![](https://user-images.githubusercontent.com/25981629/86619983-d44fa380-bfbb-11ea-91b5-d8ee4f860cf2.png)
![](https://user-images.githubusercontent.com/25981629/86620084-04974200-bfbc-11ea-8f65-11fab782a559.png)


## Results with adaptive weigths, s = 0,3,6 and r = 2 (s = sparsity thresholding parameter by which to shrink the t-statistics)
![](https://user-images.githubusercontent.com/25981629/86621615-c3ecf800-bfbe-11ea-8e9e-0b116d461b82.png)
![](https://user-images.githubusercontent.com/25981629/86621681-e2eb8a00-bfbe-11ea-8bdd-964b87325262.png)
![](https://user-images.githubusercontent.com/25981629/86621736-01518580-bfbf-11ea-80d8-c8004234aa5f.png)