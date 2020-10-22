# Stardust statistical evaluation
The figure below describe the CSSs distribution of 5 Stardust configurations (Lymph Node dataset). For the configuration '0.5' (i.e. distance based on space weights half of the distance based on transcriptional profiles) a statistical evaluation has been performed.

<img src="https://user-images.githubusercontent.com/25981629/96877431-e9bce080-1479-11eb-990a-20a0d59263a1.png" width="800" />

This script describes 100 clustering repetions. In each one the spot positions have been suffled and then Stardust has been applied.
````R
library(rCASC)
for(i in 101:200){
  system(paste("mkdir -p /media/data/users/gmotterl/LymphNode/null/",i,"/scratch",sep=""))
  setwd(paste("/media/data/users/gmotterl/LymphNode/null/",i,sep=""))
  system("cp /media/data/users/gmotterl/LymphNode/filtered_expression_matrix.txt .")
  system("cp /media/data/users/gmotterl/LymphNode/spot_coordinates.txt .")
  
  scratch.folder <- paste(getwd(),"/scratch",sep="")
  file <- paste(getwd(),"/filtered_expression_matrix.txt",sep="")
  tissuePosition <- paste(getwd(),"/spot_coordinates.txt",sep="")
  shuf <- read.delim(tissuePosition, row.names=1)
  rownames(shuf) = sample(rownames(shuf))
  write.table(shuf, tissuePosition, append = FALSE, sep = "\t",row.names = TRUE, col.names = TRUE)
  spatial2020SeuratPermutation(group="docker",scratch.folder = scratch.folder,nPerm=80,
    file=file, tissuePosition=tissuePosition, profileDistance=2, spotDistance=2, 
    spotDistanceTransformationWeight=0.5, permAtTime=8, percent=10, separator="\t",
    logTen=0, pcaDimensions=5, seed=111)
  cluster.path <- paste(data.folder=dirname(file), "Results", strsplit(basename(file),"\\.")[[1]][1], sep="/")
  cluster <- as.numeric(list.dirs(cluster.path, full.names = FALSE, recursive = FALSE))
  permAnalysisSeurat(group="docker",scratch.folder = scratch.folder,file=file, nCluster=cluster,separator="\t",sp=0.8)
  system("rm filtered_expression_matrix.txt")
  system("rm Results/filtered_expression_matrix.txt")
  system("rm Results/spot_coordinates.txt")
}

````
Then we apply Wilcoxon and Kolmogorov statistical tests.  
We compare CSSs distributions between Stardust applied to the original dataset and Stardust applied to the same dataset with shuffled positions.  
Let x be the CSSs distribution for the original dataset and y be the CSSs distribution for the k-th shuffled dataset.  

Wilcoxon test:  
H0 = true location shift between x and y is equal to 0   
H1 = true location shift between x and y is greater than 0

Kolmogorov test:   
H0 =  the CDF of x does not lie below that of y  
H1 = the CDF of x lies below that of y  (and thus x is greater than y)  


P-values are corrected with Bonferroni for multiple hypothesis testing
````R
wilcox = 1:100
kolmogorov = 1:100
zerocinque = read.delim("~/Desktop/Margarita/LymphNode/0.5/clustering_test/Results/filtered_expression_matrix/19/filtered_expression_matrix_scoreSum.txt", header=FALSE, row.names=1)
for(i in 101:200){
  path = paste("~/Desktop/Margarita/LymphNode/null/",i,"/Results/filtered_expression_matrix/",sep="")
  setwd(path)
  clusters = strsplit(list.dirs(".",recursive = FALSE)[1],"/")[[1]][2]
  path = paste(path,"/",clusters,"/filtered_expression_matrix_scoreSum.txt",sep="")
  a <- read.delim(path, header=FALSE, row.names=1)
  k = wilcox.test(zerocinque$V2, y = a$V2 ,alternative = "greater")
  wilcox[i-99] = k$p.value
  kk = ks.test(zerocinque$V2, y = a$V2 ,alternative = "less")
  kolmogorov[i-99] = kk$p.value
}
print(paste("Wilcoxon: in ", sum(wilcox < (0.05/100)), " times over 100 x is statistically higher w.r.t. to y.",sep=""))
print(paste("Kolmogorov: in ", sum(kolmogorov < (0.05/100)), " times over 100 the CDF of x lies below w.r.t. the CDF of y.",sep=""))

# Output:
# "Wilcoxon: in 100 times over 100 x is statistically higher w.r.t. to y"
# "Kolmogorov: in 100 times over 100 the CDF of x lies below w.r.t. the CDF of y."
````

In conclusion the increase in CSSs for "0.5" configuration of Stardust is statistically significative if compared with the results based on shuffled position of spots.