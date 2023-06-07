#####
library("DESeq2")
library("pheatmap")
library("data.table")
library("tibble")
edata<-fread("A://HCC/Bulkseq/CIBERSORTx/myownfiles/hcc.set.combatseq.counts.txt",header = T,check.names = F)
edata<-edata%>%column_to_rownames("Gene")
pheno<-read.table("A://HCC/Bulkseq/CIBERSORTx/myownfiles/phen.txt",header = T,check.names = F)
pheno_T<-pheno[pheno$Group=="T",]
edataT<-edata[,pheno_T$Samples]
#
classes<-read.table("A://HCC/Bulkseq/survival_class/HCC_mix_5.txt")
all.equal(pheno_T$Samples,classes$SampleID)
#
pheno_T<-cbind(pheno_T,classes$catgory)
colnames(pheno_T)[4]<-"catgory"
#
pheno_T$Group<-factor(pheno_T$Group,levels = c("N","T"))
pheno_T$Batch<-factor(pheno_T$Batch,levels = c("LOCAL","LIHC"))

###
for(i in unique(pheno_T$catgory)){
  pheno_T$catgory<-relevel(factor(pheno_T$catgory),ref = i)
  #
  ddsFullCountTable <- DESeqDataSetFromMatrix(countData = edataT,colData = pheno_T,design = ~ Batch + catgory)
  dds<-DESeq(ddsFullCountTable)
  vsd <- vst(dds,blind=F)
  save("ddsFullCountTable","dds","vsd",file = paste0("A://HCC/Bulkseq/degresult/fiveclasses/",i,".RData"))
  for(j in resultsNames(dds)[3:length(resultsNames(dds))]){
    res05 <- results(dds,name=j,alpha = 0.05)
    DEGtable<-as.data.frame(res05)
    DEGtable<-cbind(rownames(DEGtable),DEGtable)
    colnames(DEGtable)[1]<-"Gene"
    DEGtable<-DEGtable[!is.na(DEGtable$padj)&abs(DEGtable$log2FoldChange)>1&DEGtable$padj<0.05,]
    write.table(DEGtable,file = paste0("A://HCC/Bulkseq/degresult/fiveclasses/",i,"_",j,".txt"),quote = F,sep = "\t",row.names = F,col.names = T)
  }
  rm(ddsFullCountTable,dds,vsd)
}
## Fig1C - right
load("A://HCC/Bulkseq/degresult/fiveclasses/pro_T.RData")
color_select<-function(class_){
  color<-as.data.frame(readxl::read_xlsx("A://HCC/配色.xlsx",sheet = "Sheet1"))[,1:2]
  colnames(color)<-c("clusters","color")  
  loc<-match(class_,color$clusters)
  target<-color[loc,"color"]
  names(target)<-color[loc,"clusters"]
  return(target)
}
p1<-plotPCA(vsd,intgroup=c("catgory"))
score<-p1$data
p2<-ggplot(score,aes(x=PC1,y=PC2))+
  geom_point(aes(color=score[,"catgory"]),size=1)+scale_color_manual(values=color_select(class_ = score$catgory),name="catgory")+stat_ellipse(aes(fill=score[,"catgory"]),type="norm",geom="polygon",alpha=0.2,color=NA)+scale_fill_manual(values=color_select(class_ = levels(score$catgory)))+theme_classic()
#ggsave(filename = "A://HCC/figs/bulk/PCA/plotpcavsd.pdf",p2,width = 7,height = 6)

