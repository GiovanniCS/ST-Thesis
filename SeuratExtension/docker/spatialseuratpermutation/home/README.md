# Seurate extension in rCASC
This script allows to execute my extension of Seurato clustring inside rCASC container architecture.

Workflow:  
Bash
````bash
mkdit test && cd test
mkdir scratch
cp (projRootDir)/Datasets/HumanLymphNodeDocker/filtered_expression_matrix.txt .
cp (projRootDir)/Datasets/HumanLymphNodeDocker/spot_coordinates.txt .
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
spatialSeuratPermutation(group="docker",scratch.folder=paste(workingDir,"/scretch",sep=""),
    file=paste(workingDir,"/filtered_expression_matrix.txt"sep=""),
    tissuePosition=paste(workingDir,"/spot_coordinates.txt",sep=""),
    spotDistanceTransformationWeight=0.3,nPerm=2, permAtTime=2, percent=10,
    separator="\t", logTen=0, pcaDimensions=5, seed=111)

````
