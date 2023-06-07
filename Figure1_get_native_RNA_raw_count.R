library("DESeq2")
library("gplots")
library("ggplot2")
library("dichromat")
library("dplyr")
library(pheatmap)
library("org.Hs.eg.db")
library("dplyr")
library("plyr")
library(clusterProfiler)
library(fgsea)
library(data.table)
library(ggplot2)
library("reshape2")
## generate native samples raw count matrix
# source file
countloc<-"A://HCC/Bulkseq/countfiles/"
file.names.all<-sort(factor(dir(countloc,recursive = T),levels = paste0(c("CHT_LN","JKC_RN","ZFF_RN","ZMK_RN","ZYG_RN","CHT_LT","JKC_LT","ZFF_LT","ZMK_LT","ZYG_LT","CHT_RT","ZFF_RT","ZMK_RT","ZYG_RT"),".count.txt")))
sample.names.all<-unlist(lapply(strsplit(as.character(file.names.all),split = "\\.",perl = T),function(x){x[1]}))
# meta info
clinic_meta<-as.data.frame(readxl::read_xlsx("A://HCC/Seqinfo/clinic_WGCNA.xlsx"))
# construct obj
sample.names_se <- sample.names.all[sample.names.all!="ZFF_RN"]
clinic_meta_se<-clinic_meta[clinic_meta$Samples%in%sample.names_se,]
clinic_meta_se$Samples <- factor(clinic_meta_se$Samples,levels = sample.names_se)
clinic_meta_se<-arrange(clinic_meta_se,Samples)
#
conditions.T_B<-factor(ifelse(clinic_meta_se$Vascular_tumor_thrombus[match(sample.names_se,clinic_meta_se$Samples)]==1,"TBY","TBN"))
conditions_se <- conditions.T_B
file.names_se <- factor(file.names.all[file.names.all!="ZFF_RN.count.txt"])
sampleTable <- data.frame(sampleName=sample.names_se,fileName=file.names_se, condition=conditions_se)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=countloc, design=~condition)
# QC
orignID<-rownames(ddsHTSeq)
#convert
#GRCH37p13<-as.data.frame(rtracklayer::import('A://HCC/Bulkseq/gtffile/Homo_sapiens.GRCh37.p13.gtf'))
#save(GRCH37p13,file = "A://HCC/Bulkseq/gtffile/Homo_sapiens.GRCh37.p13.gtf.RData")
#
load("A://HCC/Bulkseq/gtffile/Homo_sapiens.GRCh37.p13.gtf.RData")
gtf_used<-subset(GRCH37p13,type=="gene")
gtf_used<-dplyr::select(gtf_used,gene_id,gene_name)
gtf_used$gene_id<-unlist(lapply(strsplit(x = gtf_used$gene_id,split = "\\."),function(x){x[[1]]}))
#genes which doesn't in gtf files wiil be dropped
ddsHTSeq<-ddsHTSeq[orignID%in%gtf_used$gene_id,]
map1<-mapvalues(x = rownames(ddsHTSeq),from = gtf_used$gene_id,to = gtf_used$gene_name)
##avoid duplicated genes
counts_ddsHTSeq<-as.data.frame(counts(ddsHTSeq))
counts_ddsHTSeq$gene<-map1
index=order(rowMeans(counts_ddsHTSeq[,-ncol(counts_ddsHTSeq)]),decreasing = T)
expr_ordered=counts_ddsHTSeq[index,]
keep=!duplicated(expr_ordered$gene)
expr_max=expr_ordered[keep,]
rownames(expr_max)<-expr_max$gene
ddsHTSeq_unique<-DESeqDataSetFromMatrix(countData = expr_max[,-grep("gene",colnames(expr_max))],colData = as.data.frame(colData(ddsHTSeq)),design =~condition) 
# length(unique(map1)) table(keep)
ddsHTSeq_unique <- ddsHTSeq_unique[rowSums(counts(ddsHTSeq_unique))>=10,]
dds <-DESeq(ddsHTSeq_unique)
# save counts file
merge_count<-counts(dds)
merge_count<-cbind(rownames(merge_count),merge_count)
#write.table(merge_count,"A://HCC/Bulkseq/mergecount/merge_count.txt",col.names = T,row.names = F,sep = "\t",quote = F)
merge_count <- read.table("A://HCC/Bulkseq/mergecount/merge_count.txt")
writexl::write_xlsx(x = merge_count,path = "A://HCC/submit/0606submit/elife/sourcedataforfigure/Figure1_S1/data/native_samples_raw_count.xlsx",col_names = F)
