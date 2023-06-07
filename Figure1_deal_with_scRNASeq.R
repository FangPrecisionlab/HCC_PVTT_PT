library("Seurat")
library("data.table")
library("dplyr")
library("ggplot2")
#
se_counts<-as.data.frame(fread("/xtdisk/fangxd_group/fangzhj/HCC/scRNAseq/GSE149614/GSE149614_HCC.scRNAseq.S71915.count.txt",sep = "\t",check.names = F))
rownames(se_counts)<-se_counts[,1]
se_counts<-se_counts[,-1]
metadata<-read.table("/xtdisk/fangxd_group/fangzhj/HCC/scRNAseq/GSE149614/HCC.metadata.txt",sep = "\t",header = T,row.names = 1)
#

se_obj <- CreateSeuratObject(counts = se_counts,min.cells = 3,min.features = 200,meta.data = metadata) %>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData()
#
se_obj <- RunPCA(se_obj, features = VariableFeatures(se_obj))
se_obj <- FindNeighbors(se_obj, dims = 1:30) 
se_obj <- FindClusters(se_obj, resolution = 3)
se_obj <- RunUMAP(se_obj, dims = 1:30)
saveRDS(se_obj,"/xtdisk/fangxd_group/fangzhj/HCC/scRNAseq/GSE149614/GSE149614.rds")
#
se_obj<-readRDS("A://HCC/publicdata/scRNAseq/GSE149614/GSE149614.rds")
samples<-levels(se_obj@meta.data$orig.ident)
for (i in samples) {
  tmp<-subset(se_obj,orig.ident==i)
  saveRDS(object = tmp,file = paste0("A://HCC/publicdata/scRNAseq/GSE149614/eachsample/",i,".rds"))
}
#read file
rds_loc<-"/xtdisk/fangxd_group/fangzhj/HCC/scRNAseq/GSE149614/eachsample/"
rds_file_loc <-dir(rds_loc,recursive = T,full.names = T)
rdslist<-list()
for(i in rds_file_loc){
  rdslist[[substr(i,nchar(i)-9,nchar(i)-4)]]<-readRDS(i)
}
#merge all the rds
for (i in 1:length(rdslist)) {
  rdslist[[i]] <- NormalizeData(rdslist[[i]])#
  rdslist[[i]] <- FindVariableFeatures(rdslist[[i]], selection.method = "vst")#
}
#
HCC_pub.anchors <- FindIntegrationAnchors(object.list = rdslist)##
HCC_pub <- IntegrateData(anchorset = HCC_pub.anchors)##
HCC_pub <- ScaleData(HCC_pub, features = VariableFeatures(HCC_pub))#
HCC_pub <- RunPCA(HCC_pub, features = VariableFeatures(HCC_pub))#
pc.num=1:30
HCC_pub <- FindNeighbors(HCC_pub, dims = pc.num) #
HCC_pub <- FindClusters(HCC_pub, resolution = 3)#
HCC_pub <- RunUMAP(HCC_pub, dims = pc.num)
saveRDS(HCC_pub,file = "/xtdisk/fangxd_group/fangzhj/HCC/scRNAseq/GSE149614/HCC_pub_CCA.rds")
#
HCC_pub<-readRDS("A://HCC/publicdata/scRNAseq/GSE149614/HCC_pub_CCA.rds")
HCC_pub<-RunTSNE(object = HCC_pub,dims = pc.num)
DimPlot(HCC_pub,reduction = "tsne",label.size = 10,label = T)
#saveRDS(object = HCC_pub,file = "A://HCC/publicdata/scRNAseq/GSE149614/HCC_pub_CCA.rds")
HCC_pub_markers<-FindAllMarkers(HCC_pub)
#marker gene 气泡图
library("tidyverse")

HCC_pub@meta.data$res.3<-factor(HCC_pub@meta.data$res.3,levels = c(3,4,12,15,17,19,22,24,27,29,42,43,45,47,49,1,2,8,9,11,13,14,18,20,25,28,32,35,51,5,6,10,16,21,23,26,38,39,41,44,46,52,53,34,36,37,40,50,7,30,48,31,33))
puremarker<-read.table("A://HCC/publicdata/scRNAseq/GSE149614/mainmarkers.txt",sep = "\t")
puremarker<-rev(unique(puremarker$V1))
p1<-DotPlot(HCC_pub, features = puremarker,group.by = "res.3",assay = "RNA",col.min = 0.5)+scale_x_discrete("")+scale_y_discrete("")+theme(axis.text=element_text(size=16,face="bold"))+RotatedAxis()+coord_flip()
#
Hep<-c(3,4,12,15,17,19,22,24,27,29,42,43,45,47,49)
T_NK<-c(1,2,8,9,11,13,14,18,20,25,28,32,35,51)
Mye<-c(5,6,10,16,21,23,26,38,39,41,44,46,52,53)
B_cell<-c(34,36,37,40,50)
Endo<-c(7,30,48)
Fib<-c(31,33)
current.cluster.ids<-c(Hep,T_NK,Mye,B_cell,Endo,Fib)
new.cluster.ids<-c(rep("Hep",length(Hep)),rep("T_NK",length(T_NK)),rep("Mye",length(Mye)),rep("B_cell",length(B_cell)),rep("Endo",length(Endo)),rep("Fib",length(Fib)))
HCC_pub@meta.data$celltype<-plyr::mapvalues(x = HCC_pub@meta.data$res.3,from = current.cluster.ids,to = new.cluster.ids)
DimPlot(object = HCC_pub,group.by = "celltype",label = T,label.size = 10,reduction = "tsne")
DimPlot(object = HCC_pub,group.by = "celltype",label = T,label.size = 10,reduction = "umap")
#
subtypeinfo<-read.xlsx("A://HCC/publicdata/scRNAseq/GSE149614/subtypeinfo.xlsx",sheetIndex = 1)
subtypeinfo<-subtypeinfo[,c(1,2,3,4)]
subtypeinfo<-subtypeinfo[-c(1,3),]
colnames(subtypeinfo)<-c("cell_clusters","celltype","subtype","cellnum")
subtypeinfo<-subtypeinfo[-1,]
HCC_pub@meta.data$"cell_clusters"<-paste0("C",HCC_pub@meta.data$res.3)
HCC_pub@meta.data$"subtype"<-plyr::mapvalues(x = HCC_pub@meta.data$cell_clusters,from = subtypeinfo$cell_clusters,to = subtypeinfo$subtype)
#single cell reference sample file.
set.seed(10)
HCC_pub_sub<-HCC_pub[,sample(1:ncol(HCC_pub),size = floor(round(ncol(HCC_pub))/5))]

SRS_file<-HCC_pub_sub@assays$RNA@counts
candiname<-plyr::mapvalues(x = colnames(SRS_file),from = rownames(HCC_pub_sub@meta.data),to = HCC_pub_sub@meta.data$subtype)
colnames(SRS_file)<-candiname
test<-as.data.frame(SRS_file[1:1000,1:1000])
data.table::fwrite(x = as.data.frame(SRS_file),file = "A://HCC/publicdata/scRNAseq/GSE149614/GSE149614_HCC_SRS_file.txt",sep = "\t",row.names = T,col.names = T,quote = F)
#
SRS_file<-fread("A://HCC/publicdata/scRNAseq/GSE149614/GSE149614_HCC_SRS_file.txt")
colnames(SRS_file)[1]<-"GeneSymbol"
data.table::fwrite(x = as.data.frame(SRS_file),file = "A://HCC/publicdata/scRNAseq/GSE149614/GSE149614_HCC_SRS_file.txt",sep = "\t",row.names = F,col.names = T,quote = F)
#HCC_pub<-readRDS("A://HCC/publicdata/scRNAseq/GSE149614/HCC_pub_CCA.rds")
