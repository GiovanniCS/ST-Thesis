 

 euc.dist = function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
 #//////////////////////////////////////////////////////////////////////////
 silhouette=function(nCluster,clustering.output){
 
dataPlot=cbind(as.numeric(clustering.output[,3]),as.numeric(clustering.output[,4])) 
nCluster=length(unique(clustering.output[,2]))
mainVector=as.numeric(clustering.output[,2])
intraScore=c()   
extraScore=c()
neighbor=c()
silhouetteValue=c()
    #per ogni cluster
  
    #per ogni elemento
for(k in 1:(length(dataPlot)/2))
{
    a=0
    count=0
    #per ogni altro elemento nel suo cluster
    for(j in 1:(length(dataPlot)/2)){
        if(mainVector[k]==mainVector[j])
            {   
                
                if(k != j ){
                    a=a+euc.dist(dataPlot[k,],dataPlot[j,])
                    count=count+1
                }
            }
    
    }
      intraScore[k]=a/count 



}
    extraScoreTemp=c()
extraCountTemp=c()

for(k in 1:(length(dataPlot)/2))
{
    for(s in 1:nCluster){
        extraScoreTemp[s]=0
        extraCountTemp[s]=0
    }
    
    for(j in 1:(length(dataPlot)/2))
    {
        if(mainVector[k] != mainVector[j]){
                                            extraScoreTemp[mainVector[j]]=extraScoreTemp[mainVector[j]]+ euc.dist(dataPlot[k,],dataPlot[j,])
                                            extraCountTemp[mainVector[j]]=extraCountTemp[mainVector[j]]+1
                                            }
    
    }
    extraScoreTemp=extraScoreTemp[-mainVector[k]]
    extraCountTemp=extraCountTemp[-mainVector[k]]
    extraScore[k]=min(extraScoreTemp/extraCountTemp)
    minIndex=which.min(extraScoreTemp/extraCountTemp)
    if(minIndex>=mainVector[k]){neighbor[k]=minIndex+1}else{neighbor[k]=minIndex}
    }
    for(u in 1:length(extraScore)){silhouetteValue[u]=(extraScore[u]-intraScore[u])/max(extraScore[u],intraScore[u])}
    silhouette=matrix(cbind(extraScore,intraScore,mainVector,neighbor,silhouetteValue),nrow=length(extraScore))
    colnames(silhouette) = c("extraScore","intraScore","ClusterBelong","Neighbor","SilhouetteValue") # the first row will be the header
return(cbind(clustering.output,extraScore,intraScore,neighbor,silhouetteValue))
 }
  #//////////////////////////////////////////////////////////////////////////

 recorsiveSIMLR=function(matrixCount,nCluster){
 tt=try(SIMLR(matrixCount,nCluster,cores.ratio=0))
 if(class(tt)=="try-error"){
tt=recorsiveSIMLR(matrixCount,nCluster)
			}
return(tt)
 } 
 
  recorsivegriph=function(countMatrix,image.format="pdf",use.par = FALSE){
 tt=try(griph_cluster(countMatrix,image.format="pdf",use.par = FALSE))
 if(class(tt)=="try-error"){
 
tt=recorsivegriph(countMatrix,image.format="pdf",use.par = FALSE)
			}
return(tt)
 } 
 
 
  #//////////////////////////////////////////////////////////////////////////
 #5.112772e-14

 
 
 
 
simlrF=function(matrixCount,nCluster){
cluster_result=recorsiveSIMLR(matrixCount,nCluster)
return(cbind(CellName=colnames(matrixCount),Belonging_Cluster=cluster_result$y$cluster,xChoord=cluster_result$ydata[,1],yChoord=cluster_result$ydata[,2]))
}


simlrF2=function(matrixCount,nCluster,index){
cluster_result=recorsiveSIMLR(matrixCount,nCluster)

rank=SIMLR_Feature_Ranking(cluster_result$S,matrixCount)
foo=cbind(rank$aggR,rank$pval)
foo=foo[order(foo[,1]),]
rownames(foo)=rownames(matrixCount)
foo=foo[,-1]
write.table(foo,paste("./Permutation/pvalue/",index,".",format,sep=""),sep=separator,col.names=FALSE)

return(cbind(CellName=colnames(matrixCount),Belonging_Cluster=cluster_result$y$cluster,xChoord=cluster_result$ydata[,1],yChoord=cluster_result$ydata[,2]))
}

griphF=function(countMatrix){
library("griph")
cluster_result=recorsivegriph(as.matrix(countMatrix),image.format="pdf",use.par = FALSE)
return(cbind(CellName=colnames(countMatrix),Belonging_Cluster=cluster_result$MEMB,xChoord=cluster_result$plotLVis[,1],yChoord=cluster_result$plotLVis[,2]))

}


tsneF=function(countMatrix,nCluster,perplexity){
library(Rtsne)
ts=Rtsne(dist(t(countMatrix)), is_distance=TRUE, perplexity=perplexity, verbose = TRUE)
cluster_result=kmeans(ts$Y,nCluster)
return(cbind(CellName=colnames(countMatrix),Belonging_Cluster=cluster_result$cluster,xChoord=ts$Y[,1],yChoord=ts$Y[,2]))


}




 

 
 
clustering=function(matrixName,tissuePositionFile,profileDistance,spotDistance,
    spotDistanceTransformation,nPerm,permAtTime,percent,nCluster,logTen,format,
    separator,pcaDimensions){
    if(separator=="tab"){separator2="\t"}else{separator2=separator} #BUG CORRECTION TAB PROBLEM 
    if(sparse=="FALSE"){
        countMatrix=read.table(paste("./../",matrixName,".",format,sep=""),sep=separator2,header=TRUE,row.name=1)
        if(logTen==1){countMatrix=10^(countMatrix)}
    }else{
        countMatrix <- Read10X(data.dir = "./")
        if(logTen==1){
            stop("Sparse Matrix in Seurat has to be raw count")
        }
    }
    pbmc=CreateSeuratObject(countMatrix)
    mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
    percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

    # AddMetaData adds columns to object@meta.data, and is a great place to
    # stash QC stats
    pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
    #pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
    #    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
        scale.factor = 10000)
    pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5,do.plot = FALSE)
    pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))
    pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = FALSE, pcs.print = 1:5, 
        genes.print = 5)
        
    pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
    pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE,num.pc=pcaDimensions)

DimElbowPlot2=function (object, reduction.type = "pca", dims.plot = 20, xlab = "", 
    ylab = "", title = ""){
    data.use <- GetDimReduction(object = object, reduction.type = reduction.type, 
        slot = "sdev")
    if (length(data.use) == 0) {
        stop(paste("No standard deviation info stored for", reduction.type))
    }
    if (length(x = data.use) < dims.plot) {
        warning(paste("The object only has information for", 
            length(x = data.use), "PCs."))
        dims.plot <- length(x = data.use)
    }
    data.use <- data.use[1:dims.plot]
    dims <- 1:length(x = data.use)
    data.plot <- data.frame(dims, data.use)
    return(cbind(dims,data.use))

}

PCElbowPlot2= function (object, num.pc = 20){
    return(DimElbowPlot2(object = object, reduction.type = "pca", 
        dims.plot = num.pc))
}
sdPC=PCElbowPlot2(pbmc)[,2]
#pdf("PCE_bowPlot.pdf")
#PCElbowPlot(object = pbmc)
#dev.off()

if(pcaDimensions==0){
pcaDimensions=which.max(abs(diff(sdPC)-mean(diff(sdPC))))
if(pcaDimensions==1){pcaDimensions=2}
}

pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use =seq(1,pcaDimensions), 
    resolution = 0.6, print.output = 0, save.SNN = TRUE)
    
    pbmc <- RunTSNE(object = pbmc, dims.use =seq(1,pcaDimensions), do.fast = TRUE)
    #TSNEPlot(object = pbmc)
mainVector=as.numeric(pbmc@ident) 
Coordinates=pbmc@dr$tsne@cell.embeddings

nCluster=max(mainVector)
dir.create(paste("./",nCluster,sep=""))
dir.create(paste("./",nCluster,"/Permutation",sep=""))
setwd(paste("./",nCluster,sep=""))

clustering.output= cbind(rownames(Coordinates),mainVector,Coordinates[,1],Coordinates[,2])
clustering.output=silhouette(length(unique(mainVector)),clustering.output)
colnames(clustering.output)=c("cellName","Belonging_Cluster","xChoord","yChoord","extraScore","intraScore","neighbor","silhouetteValue")
write.table(clustering.output,paste(matrixName,"_clustering.output.",format,sep=""),sep=separator2, row.names = F)
cycles=nPerm/permAtTime
cat(getwd())
for(i in 1:cycles){
    system(paste("for X in $(seq ",permAtTime,")
do
 nohup Rscript ./../../../home/permutation.R ",percent," ",matrixName," ",tissuePositionFile," ",profileDistance," ",spotDistance," ",spotDistanceTransformation," ",format," ",separator," ",logTen," ",pcaDimensions," ",sparse," $(($X +",(i-1)*permAtTime," )) & 

done"))
d=1
while(length(list.files("./Permutation",pattern=paste("*.",format,sep="")))!=i*permAtTime*2){
if(d==1){cat(paste("Cluster number ",nCluster," ",((permAtTime*i))/nPerm*100," % complete \n"))}
d=2
}

system("echo 3 > /proc/sys/vm/drop_caches")
system("sync")
gc()
}


#write.table(as.matrix(sapply(list.files("./Permutation/",pattern="cluster*"),FUN=function(x){a=read.table(paste("./Permutation/",x,sep=""),header=TRUE,col.names=1,sep=separator2)[[1]]}),col.names=1),paste(matrixName,"_",nCluster,"_clusterP.",format,sep=""),sep=separator2,row.names=FALSE, quote=FALSE)
#write.table(as.matrix(sapply(list.files("./Permutation/",pattern="killC*"),FUN=function(x){a=read.table(paste("./Permutation/",x,sep=""),header=TRUE,col.names=1,sep=separator2)[[1]]}),col.names=1),paste(matrixName,"_",nCluster,"_killedCell.",format,sep=""),sep=separator2,row.names=FALSE, quote=FALSE)

cluster_p=sapply(list.files("./Permutation/",pattern="cluster*"),FUN=function(x){a=read.table(paste("./Permutation/",x,sep=""),header=TRUE,col.names=1,sep=separator2)[[1]]})
killedC=sapply(list.files("./Permutation/",pattern="killC*"),FUN=function(x){a=read.table(paste("./Permutation/",x,sep=""),header=TRUE,col.names=1,sep=separator2)[[1]]})

write.table(as.matrix(cluster_p,col.names=1),paste(matrixName,"_",nCluster,"_clusterP.",format,sep=""),sep=separator2,row.names=FALSE, quote=FALSE)
write.table(as.matrix(killedC,col.names=1),paste(matrixName,"_",nCluster,"_killedCell.",format,sep=""),sep=separator2,row.names=FALSE, quote=FALSE)


pdf("hist.pdf")
clusters=apply(cluster_p,2,FUN=function(x){max(x)})
hist(clusters,xlab="nCluster",breaks=length(unique(cluster_p)))
dev.off()

write.table(sort(unique(clusters)),paste("./../rangeVector.",format,sep=""),sep=separator2,row.names=FALSE,col.names=FALSE)
system("rm -r Permutation")
return(length(unique(mainVector)))

}
relationMatrix=function(mainVector,nameVector){
 rel.matrix=as.numeric(mainVector) %*% t(1/as.numeric(mainVector))
 rel.matrix[rel.matrix!=1]=0
 colnames(rel.matrix)=nameVector
 rownames(rel.matrix)=nameVector
 return(rel.matrix)
 
}



silhouettePlot=function(matrixName,rangeVector,format,separator){
if(separator=="tab"){separator="\t"} #BUG CORRECTION TAB PROBLEM 
count=1
l=list()
for(i in rangeVector){
l[[count]]=read.table(paste("./",i,"/",matrixName,"_clustering.output.",format,sep=""),sep=separator,header=TRUE)[,8]
count=count+1
}
pdf(paste(matrixName,"_vioplot.pdf",sep=""))
do.call(vioplot,c(l,list(names=rangeVector)))
dev.off()
}
