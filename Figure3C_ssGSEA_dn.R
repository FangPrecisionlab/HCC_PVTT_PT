#Make input files log TPM 
Gene_Exp<-read.table("A://HCC/Bulkseq/WGCNA/HCC_mix/mix_TPM_T_log2.txt")
GCT_file<-rbind(NA,NA,Gene_Exp)
GCT_file<-cbind(NA,GCT_file)
GCT_file[1,1]<-"#1.2"
GCT_file[2,c(1,2)]<-c(nrow(Gene_Exp)-1,ncol(Gene_Exp)-1)
GCT_file[3:nrow(GCT_file),1]<-GCT_file[3:nrow(GCT_file),2]
GCT_file[3,c(1,2)]<-c("NAME","Description")
write.table(GCT_file,"A://HCC/Bulkseq/ssGSEA/5classes/TCGA_LIHC_T.gct",sep = "\t",quote = F,na = "",col.names = F,row.names = F)
#go to ssgsea 
#https://cloud.genepattern.org/gp/pages/index.jsf?lsid=urn:lsid:broad.mit.edu:cancer.software.genepattern.module.analysis:00270:10.1.0
#
scores<-read.table("A://HCC/Bulkseq/ssGSEA/5classes/stemness_5_classes.gct",sep = "\t",fill = T)
scores<-read.table("A://HCC/Bulkseq/ssGSEA/5classes/5classes_cell_proli.gct",sep = "\t",fill = T)
scores<-read.table("A://HCC/Bulkseq/ssGSEA/5classes/5classes_metas.gct",sep = "\t",fill = T)
scores<-read.table("A://HCC/Bulkseq/ssGSEA/5classes/exhausted_THCC.gct",sep = "\t",fill = T)
scores<-read.table("A://HCC/Bulkseq/ssGSEA/5classes/normal_liver.gct",sep = "\t",fill = T)

getbox<-function(note){
  scores<-t(scores[2:3,3:ncol(scores)])
  colnames(scores)<-c("samples","scores")
  scores<-as.data.frame(scores)
  scores$scores<-scale(as.numeric(scores$scores))
  #
  meta<-read.table("A://HCC/Bulkseq/survival_class/HCC_mix_5.txt")
  Metas<-rownames(meta)[which(meta$catgory=="pro_meta")]
  pro_T<-rownames(meta)[which(meta$catgory=="pro_T")]
  Mix<-rownames(meta)[which(meta$catgory=="T_meta_mix")]
  NK<-rownames(meta)[which(meta$catgory=="NKC_")]
  MemT<-rownames(meta)[which(meta$catgory=="MemT_")]
  #
  scores$"class"<-plyr::mapvalues(x = scores$samples,from = c(Metas,pro_T,Mix,NK,MemT),to = c(rep("Metas",length(Metas)),rep("pro_T",length(pro_T)),rep("Mix",length(Mix)),rep("NK",length(NK)),rep("MemT",length(MemT))))
  scores$class<-factor(scores$class,levels = c("pro_T","Mix","Metas","NK","MemT"))
  #
  t.test(scores[scores$class=="Metas","scores"],scores[scores$class=="pro_T","scores"])
  #
  library("ggplot2")
  library(ggpubr)
  my_comparisons <- list(c("pro_T","Mix"),c("Mix","Metas"),c("pro_T","Metas"))
  scores<-subset(scores,!class%in%c("NK","MemT"))
  p<<-ggboxplot(data = scores,x = "class",y = "scores", color = "class", title = note,
               palette = c("#00AFBB", "#E7B800", "#FC4E07","red","blue"), 
               add = "jitter", shape="class",bxp.errorbar = T,bxp.errorbar.width = 0.05,width=0.7,font.label = list(size = 17, color = "black"),label = NULL)+stat_compare_means(comparisons = my_comparisons, method = "t.test",label = c("p.signif"),size=5,paired = F,method.args = list(alternative="two.sided",paired=F))+theme(text=element_text(size=25))
  #ggsave(filename = paste0("A://HCC/figs/bulk/GSEA/ssGSEA/",note,".pdf"),plot = p,height = 10,width = 10)
}
getbox(note = "stemness")
#getbox(note = "prolif")
#getbox(note = "metas")
#getbox(note = "exhausted")
#getbox(note = "normal")
