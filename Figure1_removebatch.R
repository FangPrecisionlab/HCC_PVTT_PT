library("dplyr")
library("plyr")
library("data.table")
library("tibble")
library("convert")
library("ggplot2")
library("DESeq2")
####import data 
TCGA_counts<-fread("A://HCC/Bulkseq/TCGA_LIHC/expr_ma/rawcount/rawcount.txt")
Local_counts<-read.table("A://HCC/Bulkseq/mergecount/merge_count.txt",header = T)
##meta info 
#local
clinic_meta<-as.data.frame(readxl::read_xlsx("A://HCC/Seqinfo/clinic_WGCNA.xlsx"))
clinic_meta<-clinic_meta %>% dplyr::select("Samples","class")
clinic_meta$class<-mapvalues(x = clinic_meta$class,from = c(0,1),to = c("N","T"))
clinic_meta$"Batch"<-"LOCAL"
#tcga
Normal<-colnames(TCGA_counts)[-1][unlist(lapply(strsplit(colnames(TCGA_counts)[-1],split = "-"),function(x){substr(x[4],1,2)}%in%c(10:29)))]
Tumor<-colnames(TCGA_counts)[-1][!unlist(lapply(strsplit(colnames(TCGA_counts)[-1],split = "-"),function(x){substr(x[4],1,2)}%in%c(10:29)))]
LIHC_meta<-as.data.frame(cbind(Normal,"N")%>%rbind(cbind(Tumor,"T")))
LIHC_meta$"Batch"<-"LIHC"
colnames(LIHC_meta)<-colnames(clinic_meta)
#merge
phen<-rbind(clinic_meta,LIHC_meta)
colnames(phen)[2]<-"Group"
phen<-phen[phen$Samples!="ZFF_RN",]
rm(clinic_meta,LIHC_meta)
#write.table(x = phen,file = "A://HCC/Bulkseq/CIBERSORTx/myownfiles/phen.txt",col.names = T,row.names = F,sep="\t",quote = F)
##ExpressionSet
get_expr_Set<-function(x=phen,y=countslist,ncount=10,occurrence=0.2){
  #x=phen
  #y=countslist
  #ncount=10
  #occurrence=0.2
  prof = y[[1]]#第一个数据集
  for(i in 2:length(y)){
    prof <- inner_join(prof,y[[i]], by="gene")#保留共同基因名的counts数据
  }
  prf <- prof %>% column_to_rownames("gene")
  
  sid <- intersect(x$Samples,colnames(prf))#和counts信息适配
  phe <- x %>% filter(Samples %in% sid) %>% 
    mutate(Batch = factor(as.character(Batch)),
           Group=factor(as.character(Group))) %>% 
    mutate(batch_number = as.numeric(Batch)) %>%
    column_to_rownames("Samples")
  #QC
  prf.cln <- prf %>% rownames_to_column("tmp") %>% 
    filter(apply(dplyr::select(.,-one_of("tmp")),1,function(x){
      sum(x!=0)/length(x)}) > occurrence) %>%
    dplyr::select(c(tmp,rownames(phe)))%>%
    column_to_rownames("tmp")
  
  temp <-apply(prf.cln,1,function(x){length(unique(x[x>0]))}) %>% 
    data.frame() %>% setNames("number")
  remain_genes <- temp %>% filter(number >2)
  
  prf.cln <- prf.cln[rowSums(prf.cln)>ncount,]
  prf.cln2 <- prf.cln[rownames(prf.cln)%in%rownames(remain_genes),]
  
  #确定表达矩阵和表型的顺序
  for(i in 1:ncol(prf.cln2)){
    if(!(colnames(prf.cln2))[i]==rownames(phe)[i]){
      stop(paste0(i,"Wrong"))
    }
  }
  
  exprs <- as.matrix(prf.cln2)
  adf <- new("AnnotatedDataFrame",data=phe)
  experimentData <- new("MIAME",
                        name="Fangzj",
                        lab ="Fang",
                        contact="fangzj@big.ac.cn",
                        title="Experiment",
                        abstract="The gene ExpressionSet count",
                        other=list(note="Created from count files"))
  expressionSet<-new("ExpressionSet",
                     exprs=exprs,
                     phenoData=adf,
                     experimentData=experimentData)
  return(expressionSet)
}
pca_fun<-function(expers_set = hcc.set,color="Batch",shape="Group"){
  #expers_set = hcc.set
  pheno<- pData(expers_set)
  edata<-exprs(expers_set)
  # scale the data
  edata<-vst(edata)
  pca<-prcomp(t(edata),scale.=F,center=F)
  score <- inner_join(pca$x %>% data.frame() %>%
                        rownames_to_column("SampleID") %>%
                        select(c(1:3)),
                      pheno %>% rownames_to_column("SampleID"),
                      by="SampleID")
  p1<-ggplot(score,aes(x=PC1,y=PC2))+
    geom_point(aes(color=score[,color],shape=score[,shape]),size=3.5)+theme_bw()+stat_ellipse(aes(fill=Group),type="norm",geom="polygon",alpha=0.2,color=NA)
  return(p1)
}
#
countslist<-list(Local_counts,TCGA_counts)
hcc.set<-get_expr_Set(y = countslist)

# remove batch
library(sva)
ComBat_seq_fun<-function(x=hcc.set){
  #x=hcc.set
  qcMetadata<-pData(x)
  qcData <- t(exprs(x)) %>% data.frame(.,check.names = F)
  #确定表达矩阵和表型的顺序
  for(i in 1:nrow(qcData)){
    if(!(rownames(qcData))[i]==rownames(qcMetadata)[i]){
      stop(paste0(i,"Wrong"))
    }
  }
  #
  adjusted_counts<-ComBat_seq(counts = t(qcData),
                              batch=qcMetadata$Batch,
                              group=qcMetadata$Group,
                              covar_mod = NULL,
                              full_mod = T,
                              shrink = F,
                              shrink.disp = F,
                              gene.subset.n = NULL)
  exprs <- as.matrix(adjusted_counts,check.names=FALSE)
  adf <- new("AnnotatedDataFrame",data=qcMetadata)
  experimentData <- new("MIAME",
                        name="Fangzj",lab ="Fang",contact="fangzj@big.ac.cn",
                        title="Experiment",
                        abstract="The gene ExpressionSet",
                        other=list(note="adjusted counts by ComBat_seq"))
  expressionSet<-new("ExpressionSet",exprs=exprs,
                     phenoData=adf,
                     experimentData=experimentData)
  return(expressionSet)
  
}
# pca 
pca_fun<-function(expers_set = hcc.set,color="Batch",shape="Group"){
  #expers_set = hcc.set
  pheno<- pData(expers_set)
  edata<-exprs(expers_set)
  # scale the data
  edata<-vst(edata)
  pca<-prcomp(t(edata),scale.=F,center=F)
  score <- inner_join(pca$x %>% data.frame() %>%
                        rownames_to_column("SampleID") %>%
                        select(c(1:3)),
                      pheno %>% rownames_to_column("SampleID"),
                      by="SampleID")
  p1<-ggplot(score,aes(x=PC1,y=PC2))+
    geom_point(aes(color=score[,color],shape=score[,shape]),size=3.5)+theme_bw()+stat_ellipse(aes(fill=Group),type="norm",geom="polygon",alpha=0.2,color=NA)
  return(p1)
}
#校正结果
hcc.set.combatseq<-ComBat_seq_fun(x = hcc.set)
cowplot::plot_grid(pca_fun(expers_set=hcc.set,color = "Group",shape = "Batch"),
                   pca_fun(expers_set=hcc.set.combatseq,color = "Group",shape = "Batch"),
                   ncol=2,
                   labels = c("origin","combatseq"))
for_ciber_count<-hcc.set.combatseq@assayData$exprs
for_ciber_count<-cbind(rownames(for_ciber_count),for_ciber_count)
colnames(for_ciber_count)[1]<-"Gene"
#write.table(x = for_ciber_count,file = "A://HCC/Bulkseq/CIBERSORTx/myownfiles/hcc.set.combatseq.counts.txt",sep = "\t",quote = F,row.names = F)
#perform CIBERSORTx analysis on for_ciber_count. Please refer to website: https://cibersortx.stanford.edu/index.php
for_ciber_count <- data.table::fread("A://HCC/Bulkseq/CIBERSORTx/myownfiles/hcc.set.combatseq.counts.txt")
writexl::write_xlsx(x = for_ciber_count,path = "A://HCC/submit/0606submit/elife/sourcedataforfigure/Figure1_S1/data/batch_removed_raw_count.xlsx")
