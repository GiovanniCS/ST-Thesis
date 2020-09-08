# rCASC clustering stability evaluation
In this experiment, the Seurat clustering extension is evaluated trough rCASC 
stability score. In my extension, the distance matrix used to compute the SNN network
(Shared Nearest Neighbor) is computed manually adding two distinct distance measures.  
The first one is the euclidean distance (default parameter "profileDistance") between
transcriptional profiles in PCA space.  
The second one is the euclidean distance (default parameter "spotDistance") between
physical distance of spots on the 10X chip.  
Moreover "spotDistanceTransformationWeight" parameter can be decided in [0,1] in order
to give a weigth to the second measure during linear scaling step. The second measure is 
scaled according to the following formula

<img src="https://user-images.githubusercontent.com/25981629/92521037-d84aac80-f21c-11ea-9d50-e43189e8d29f.png" width="800" />

Workflow:

````bash
mkdir -p clustering_test/scratch
cd clustering_test
wget https://github.com/GiovanniCS/ST-Thesis/raw/master/Datasets/HumanLymphNodeDocker/filtered_expression_matrix.txt.zip
wget https://raw.githubusercontent.com/GiovanniCS/ST-Thesis/master/Datasets/HumanLymphNodeDocker/spot_coordinates.txt
unzip filtered_expression_matrix.txt.zip
rm -rf __MACOSX/ filtered_expression_matrix.txt.zip

# skip download of unnecessary containers through manual pull
docker pull giovannics/spatial2020seuratpermutation
docker pull repbioinfo/seuratanalysis
````
````R
install.packages("devtools")
install_github("giovannics/rcasc")
library(devtools)
library(rCASC)

scratch.folder <- paste(getwd(),"/scratch",sep="")
file <- paste(getwd(),"/filtered_expression_matrix.txt",sep="")
tissuePosition <- paste(getwd(),"/spot_coordinates.txt",sep="")
spatial2020SeuratPermutation(group="docker",scratch.folder = scratch.folder,nPerm=80, 
    file=file, tissuePosition=tissuePosition, profileDistance=2, spotDistance=2, 
    spotDistanceTransformationWeight=1, permAtTime=8, percent=10, separator="\t",
    logTen=0, pcaDimensions=5, seed=111)

cluster.path <- paste(data.folder=dirname(file), "Results", strsplit(basename(file),"\\.")[[1]][1], sep="/")
cluster <- as.numeric(list.dirs(cluster.path, full.names = FALSE, recursive = FALSE))

permAnalysisSeurat(group="docker",scratch.folder = scratch.folder,file=file, nCluster=cluster,separator="\t",sp=0.8)

# Repeat with spotDistanceTransformationWeight = 0, 0.25, 0.5 and 0.75
````

Clustering visualization and Stability_Violin_Plot, spotDistanceTransformationWeight = 1 (max of position distance is equal to the max of transcriptional profiles distance)
<p float="left">
  <img src="https://user-images.githubusercontent.com/25981629/92493331-ac1b3580-f1f4-11ea-8296-c954082e2035.png" width="400" />
  <img src="https://user-images.githubusercontent.com/25981629/92493708-1df37f00-f1f5-11ea-9c11-2369fd0b95a4.png" width="300" /> 
</p>

Clustering visualization and Stability_Violin_Plot, spotDistanceTransformationWeight = 0.75 (max of position distance is 75% of the max of transcriptional profiles distance)
<p float="left">
  <img src="https://user-images.githubusercontent.com/25981629/92522204-b6eac000-f21e-11ea-94ef-52fcd0bc6357.png" width="400" />
  <img src="https://user-images.githubusercontent.com/25981629/92522022-7428e800-f21e-11ea-8413-333541de326a.png" width="300" /> 
</p>


Clustering visualization and Stability_Violin_Plot, spotDistanceTransformationWeight = 0.5 (max of position distance is 50% of max of transcriptional profiles distance)
<p float="left">
  <img src="https://user-images.githubusercontent.com/25981629/92503253-46817600-f201-11ea-92c9-7751ad5a66ff.png" width="400" />
  <img src="https://user-images.githubusercontent.com/25981629/92503048-fa363600-f200-11ea-8647-ef7ea2fa24f0.png" width="300" /> 
</p>

Clustering visualization and Stability_Violin_Plot, spotDistanceTransformationWeight = 0.25 (max of position distance is 25% of the max of transcriptional profiles distance)
<p float="left">
  <img src="https://user-images.githubusercontent.com/25981629/92512067-86028f00-f20e-11ea-863c-cdc1a046ffc0.png" width="400" />
  <img src="https://user-images.githubusercontent.com/25981629/92511889-40de5d00-f20e-11ea-9a47-1005dbf183e4.png" width="300" /> 
</p>

Clustering visualization and Stability_Violin_Plot, spotDistanceTransformationWeight = 0 ( no position data)
<p float="left">
  <img src="https://user-images.githubusercontent.com/25981629/92497868-2c906500-f1fa-11ea-8d8a-5c52a49a0fa0.png" width="400" />
  <img src="https://user-images.githubusercontent.com/25981629/92497683-f2bf5e80-f1f9-11ea-8db6-e585f67e1391.png" width="300" /> 
</p>
