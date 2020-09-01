# Seurate extension in rCASC
This script allows to execute my extension of Seurat clustring inside rCASC container architecture.

Workflow:  
Bash
````bash
mkdit test && cd test
mkdir scratch
cp (projRootDir)/Datasets/HumanLymphNodeDocker/filtered_expression_matrix.txt .
cp (projRootDir)/Datasets/HumanLymphNodeDocker/spot_coordinates.txt .
cp (projRootDir)/Datasets/HumanLymphNodeDocker/spot_coordinates_no_header.txt .
R
````
R
````R
install.packages("devtools")
library(devtools)
install_github("giovannics/rCASC", ref="master")
library(rCASC)

downloadContainers()
workingDir=getwd()
scratch.folder=paste(workingDir,"/scretch",sep="")
file=paste(workingDir,"/filtered_expression_matrix.txt"sep="")
tissuePosition=paste(workingDir,"/spot_coordinates.txt",sep="")

spatial2020SeuratPermutation(group="docker",scratch.folder=scratch.folder,
    file=file, tissuePosition=tissuePosition, spotDistanceTransformationWeight=0.3,
    nPerm=2, permAtTime=2, percent=10, separator="\t", logTen=0, pcaDimensions=5, 
    seed=111)

cluster.path <- paste(data.folder=dirname(file), "Results", 
    strsplit(basename(file),"\\.")[[1]][1], sep="/")
cluster <- as.numeric(list.dirs(cluster.path, full.names = FALSE, recursive = FALSE))

permAnalysisSeurat(group="docker", scratch.folder=scratch.folder,file=file, 
    nCluster=cluster, separator="\t", sp=0.8)

spatialAnalysis2(group="docker", scratch.folder=scratch.folder,file=file, 
    nCluster=cluster, separator="\t", tissuePosition=tissuePosition)

````