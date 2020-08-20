#!/usr/bin/env Rscript
library(Seurat)
library(RColorBrewer)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
profilesDistance = 2
spotDistance = 2
spotDistanceTransformation = "none"
workingDir = ""
print("Parsing the provided parameters ..")
for(i in args){
    param = strsplit(i,"=")[[1]]
    if(param[1] == "profilesDistance"){
        profilesDistance = suppressWarnings(as.integer(param[2]))
        if(is.na(profilesDistance))
            stop(paste("Provided param profilesDistance is not a number. ",
                paste(param,collapse="=")))
    }
    else if (param[1] == "spotDistance") {
        spotDistance = suppressWarnings(as.integer(param[2]))
        if(is.na(spotDistance))
            stop(paste("Provided param spotDistance is not a number. ",
                paste(param,collapse="=")))
    }
    else if (param[1] == "spotDistanceTransformation") {
        spotDistanceTransformation = param[2]
        if(!(spotDistanceTransformation == "linear" | spotDistanceTransformation == "log" |
            spotDistanceTransformation == "sqrt" | spotDistanceTransformation == "none"))
            stop(paste("Param spotDistanceTransformation must have one of the",
                " following values: 'none', 'linear', 'log' or 'sqrt'"))
    }
    else if (param[1] == "workingDir") {
        workingDir = param[2]
        if (!dir.exists(workingDir))
            stop(paste("Directory '",workingDir,"' not found"))
        setwd(workingDir)
    }
    else {
        stop(paste("parameter not recongnized: ",param[1]))
    }
}
if(workingDir == ""){
    stop(paste("Working directory not provided, use workingDir=path_to_working_dir as parameter"))
}
print("Loading the dataset..")
mouse <- Load10X_Spatial(data.dir = "./filtered_feature_bc_matrix")
print("Normalization of transcriptional data..")
mouse <- suppressWarnings(SCTransform(mouse, assay = "Spatial", verbose = FALSE))
print("Compiuting distance measures..")
mouse <- RunPCA(mouse, assay = "SCT", verbose = FALSE)
m = mouse@reductions[["pca"]]@cell.embeddings
distPCA = dist(m,method="minkowski",p=profilesDistance)  
coord <- GetTissueCoordinates(mouse)
distcoord = dist(coord,method="minkowski",p=spotDistance)
if (spotDistanceTransformation == "linear" & max(distcoord) > max(distPCA)){
    distcoord = distcoord*(max(distPCA)/max(distcoord))
} else if (spotDistanceTransformation == "log"){
    distcoord = log(distcoord+1)
} else if (spotDistanceTransformation == "sqrt") {
   distcoord = sqrt(distcoord)
}
finaldist = distPCA + distcoord
print("Building SSN and appling community detection..")
neighbors = FindNeighbors(finaldist)
neighbors = list(neighbors_nn=neighbors$nn,neighbors_snn=neighbors$snn)
mouse@graphs = neighbors
mouse <- suppressWarnings(FindClusters(mouse, verbose = FALSE,
    graph.name = "neighbors_snn"))
print("Generating output figure..")
palette = colorRampPalette(brewer.pal(11,name="BrBG"))
clusters = length(levels(mouse@active.ident))
jpeg(file=paste("profilesDistance:",profilesDistance," spotDistance:",
    spotDistance," spotDistanceTransformation:",spotDistanceTransformation,".jpeg"))
title.style = element_text(size=12)
title = paste("Minkowski dist. of param: ",profilesDistance," between profiles \nMinkowski dist. of param: ",
              spotDistance," between fisical spots \nFisical distance transformation: ",spotDistanceTransformation)
suppressWarnings(SpatialDimPlot(mouse, label = TRUE, label.size = 3,
    cols=palette(clusters)) + ggtitle(title) + theme(plot.title=title.style))
out = dev.off()
