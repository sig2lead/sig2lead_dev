###############################
# ChemmineR Clustering of Analog
##############################
library(ChemmineR)
#library(help="ChemmineR")
library("gplots")
library("cluster")
library("ggplot2")

cluster_compounds <- function()
  {
sdfset <- sdf_smiles
#sdfset <- read.SDFset("sdf_smiles.sdf")
#sdfset <- read.SDFset("lib/A1TestedWithiLINCSExtract")
#sdfset <- sdf_test
valid <- validSDF(sdfset)
unique_ids <- makeUnique(sdfid(sdfset))
cid(sdfset) <- unique_ids
rpropma <- data.frame(MF=MF(sdfset), MW=MW(sdfset))
#plot(sdfset[1:3], print=FALSE)


apset <-sdf2ap(sdfset)
fpset<<-desc2fp(apset)
print("Calculating Similarity Matrix")
#begin_sim_time <- Sys.time()
simMA <<- sapply(cid(fpset), function(x) fpSim(fpset[x], fpset, sorted=FALSE))

hc <<-hclust(as.dist(1-simMA), method="average")
#hc <<-hclust(as.dist(1-simMA), method="complete")
par(cex=0.3)
#plot(as.dendrogram(hc),edgePar=list(col=4, lwd=2),horiz=TRUE)

#heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5)

#cut.4<-cutree(hc, h=0.4)
#cut.2<-cutree(hc, h=0.2)

}

bypass_clustering <- function(lsm_sdf, added_sdf)
{
  #sdfset <- sdf_smiles_2
  #sdfset <- read.SDFset("sdf_smiles.sdf")
  #sdfset <- read.SDFset("lib/A1TestedWithiLINCSExtract")
  #sdfset <- sdf_test
  #valid <- validSDF(sdfset)
 # unique_ids <- makeUnique(sdfid(sdfset))
  #cid(sdfset) <- unique_ids
  #rpropma <- data.frame(MF=MF(sdfset), MW=MW(sdfset))
  #plot(sdfset[1:3], print=FALSE)
  
  
  #apset <-sdf2ap(sdfset)
  #fpset<<-desc2fp(apset)
  
  apset_lsm <-sdf2ap(lsm_sdf)
  fpset_lsm <<-desc2fp(apset_lsm)
  
  apset_added <-sdf2ap(added_sdf)
  fpset_added <<-desc2fp(apset_added)
  
  print("Starting to Calculate Similarity Matrix")
  begin_sim_calc <- Sys.time()
   #simMA <<- sapply(cid(fpset), function(x) fpSim(fpset[x], fpset, sorted=FALSE))
  simMA <<- sapply(cid(fpset_added), function(x) fpSim(fpset_added[x], fpset_lsm, sorted=FALSE))
  colnames(simMA) <- sdfid(subset_add)
  simMA <<- simMA
  #simMA <<- fpSim(x=lsm_compounds, y=added_compounds, method="Tanimoto")
  end_sim_calc <- Sys.time()
  total_sim_time <<- end_sim_calc - begin_sim_calc
  print(paste("Total Time to Calculate Similarity Matrix: ", total_sim_time, sep= " "))
  #hc <<-hclust(as.dist(1-simMA), method="average")
  #hc <<-hclust(as.dist(1-simMA), method="complete")
  par(cex=0.3)
  #plot(as.dendrogram(hc),edgePar=list(col=4, lwd=2),horiz=TRUE)
  
  #heatmap.2(1-simMA, Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc), col=colorpanel(40, "white","yellow","red"), density.info="none", trace="none", labCol=cid(fpset), labRow=cid(fpset), cexRow=0.5, cexCol=0.5)
  
  #cut.4<-cutree(hc, h=0.4)
  #cut.2<-cutree(hc, h=0.2)
  
}
####################################
#Get Clusters
####################################

#rect.hclust(hc, h=.4, border="red")

#cut.4[hc$order]
#ClusterMembers<-cbind(clusterID=cut.4)
#cluster1 <- 1-simMA[cut.4==49]
#cluster1

#ClusterMembers[hc$order,]

      