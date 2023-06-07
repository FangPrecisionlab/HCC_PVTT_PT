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


load("A://HCC/Bulkseq/degresult/fiveclasses/pro_T.RData")
# pro_meta_vs_proT padj<1e-25,abs>1
res05 <- results(dds,name=resultsNames(dds)[5],alpha = 0.05)
sig.dat <- assay(vsd)[res05$padj < 1e-25 &!is.na(res05$padj)&(abs(res05$log2FoldChange)>1),]
dim(sig.dat)
#
df <- as.data.frame(colData(dds)[,c("Batch","catgory")])
#meta_vs_proT
sig.dat_meta_vs_T<-sig.dat[,pheno_T[pheno_T$catgory%in%c("pro_T","pro_meta"),"Samples"]]
#
color_select<-function(class_){
  color<-as.data.frame(readxl::read_xlsx("A://HCC/配色.xlsx",sheet = "Sheet1"))[,1:2]
  colnames(color)<-c("clusters","color")  
  loc<-match(class_,color$clusters)
  target<-color[loc,"color"]
  names(target)<-color[loc,"clusters"]
  return(target)
}
labs.col <- colnames(sig.dat_meta_vs_T)
labs.col[grep("TCGA",labs.col)] <- ""

p2<-pheatmap::pheatmap(sig.dat_meta_vs_T,scale="row", fontsize_row=9,annotation_col = df,show_rownames = F,treeheight_row = 0,color = colorRampPalette(c("blue","white","red"))(23),show_colnames = T,annotation_colors = list(catgory=color_select(class_ = levels(df$catgory))),labels_col = labs.col)
#ggsave(filename = "A://HCC/figs/bulk/heatmap/DEGs/meta_VS_proT.pdf",plot = p2,width = 7,height = 5)