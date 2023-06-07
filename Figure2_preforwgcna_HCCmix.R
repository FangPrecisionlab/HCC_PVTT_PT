library("data.table")
library("dplyr")
library("tibble")
library("plyr")
#adjusted counts
hcc_mix_counts<-as.data.frame(fread("A://HCC/Bulkseq/CIBERSORTx/myownfiles/hcc.set.combatseq.counts.txt"))
hcc_mix_counts<-hcc_mix_counts %>% column_to_rownames("Gene")
###########Cal FPKM
library(GenomicFeatures)
id_length<-read.table("A://HCC/Bulkseq/id_gene_lens_grch37.txt",header = T)
load("A://HCC/Bulkseq/gtffile/Homo_sapiens.GRCh37.p13.gtf.RData")
gtf_used<-subset(GRCH37p13,type=="gene")
gtf_used<-dplyr::select(gtf_used,gene_id,gene_name)
gtf_used$gene_id<-unlist(lapply(strsplit(x = gtf_used$gene_id,split = "\\."),function(x){x[[1]]}))
#
used_id<-gtf_used[match(rownames(hcc_mix_counts),gtf_used$gene_name),"gene_id"]
id_length_fil <- id_length[match(used_id,id_length$gene_id),]
id_length_fil$"gene_name"<-gtf_used[match(id_length_fil$gene_id,gtf_used$gene_id),"gene_name"]
#
counts<-hcc_mix_counts
counts$"gene_name"<-rownames(counts)
mergecount <- merge(counts, id_length_fil, by = 'gene_name')
rownames(mergecount)<-mergecount$gene_name
mergecount<-mergecount[rownames(counts),]
##Cal TPM
corredcount<-mergecount[,2:(ncol(mergecount)-2)]/(mergecount[,ncol(mergecount)]/1000)
all.equal(counts$gene_name,mergecount$gene_name)
#Rtotal
Rtotal<-colSums(corredcount)
#
TPM<-t(apply(X = mergecount[,-grep("gene_name|gene_id",colnames(mergecount))],MARGIN = 1,FUN = function(x){x[-length(x)]*1000*1000000/(x[length(x)]*Rtotal)}))
#
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}
write.table(adjustdata(data = TPM),"A://HCC/Bulkseq/WGCNA/HCC_mix/HCC_mix_TPM.txt",quote = F,col.names = T,row.names = F,sep = "\t")
fwrite(x = adjustdata(data = TPM),file = "A://HCC/Bulkseq/WGCNA/HCC_mix/HCC_mix_TPM.txt",quote = F,col.names = T,row.names = F,sep = "\t")
##

TPM<-fread("A://HCC/Bulkseq/WGCNA/HCC_mix/HCC_mix_TPM.txt",header = T,check.names = F,sep = "\t")
TPM<-TPM%>%column_to_rownames("V1")
TPM_log2<-log(x = TPM+1,base = 2)
TPM_log2<-cbind(rownames(TPM_log2),TPM_log2)
colnames(TPM_log2)[1]<-"gene"
#
Normal<-colnames(TPM_log2)[-(1:14)][unlist(lapply(strsplit(colnames(TPM_log2)[-(1:14)],split = "-"),function(x){substr(x[4],1,2)}%in%c(10:29)))]
Tumor<-colnames(TPM_log2)[-(1:14)][!unlist(lapply(strsplit(colnames(TPM_log2)[-(1:14)],split = "-"),function(x){substr(x[4],1,2)}%in%c(10:29)))]
#
TPM_N_log2<-TPM_log2[,c(grep("gene|LN|RN",colnames(TPM_log2)),grep(paste0(Normal,collapse = "|"),colnames(TPM_log2)))]
TPM_T_log2<-TPM_log2[,c(grep("gene|LT|RT",colnames(TPM_log2)),grep(paste0(Tumor,collapse = "|"),colnames(TPM_log2)))]
#
write.table(TPM_N_log2,"A://HCC/Bulkseq/WGCNA/HCC_mix/mix_TPM_N_log2.txt",col.names = T,row.names = F,sep = "\t",quote = F)
write.table(TPM_T_log2,"A://HCC/Bulkseq/WGCNA/HCC_mix/mix_TPM_T_log2.txt",col.names = T,row.names = F,sep = "\t",quote = F)
write.table(TPM_log2,"A://HCC/Bulkseq/WGCNA/HCC_mix/mix_TPM_log2.txt",col.names = T,row.names = F,sep = "\t",quote = F)
##FPKM
FPKM<-as.data.frame(t(t(corredcount)/colSums(counts[-grep("gene_name",colnames(counts))]) * 10^6))
write.table(adjustdata(data = FPKM),"A://HCC/Bulkseq/WGCNA/HCC_mix/allthefpkm.txt",quote = F,col.names = T,row.names = F,sep = "\t")

####表型数据准备
library(tidyr)
library("data.table")
TPM_log<-fread("A://HCC/Bulkseq/WGCNA/HCC_mix/mix_TPM_T_log2.txt",header = T,check.names = F,sep = "\t")
TPM_log<-TPM_log%>%column_to_rownames("gene")
#
trait.Local = xlsx::read.xlsx("A://HCC/Seqinfo/clinic_WGCNA.xlsx",sheetIndex = 1,header = T)
trait.LIHC <- as.data.frame(readr::read_tsv("A://HCC/Bulkseq/TCGA_LIHC/clinical/clinical.tsv"))
#
phenotype<-read.table("A://HCC/Bulkseq/CIBERSORTx/myownfiles/phen.txt",header = T)
phenotype<-left_join(x = phenotype,y = trait.Local,by="Samples")
#
trait.LIHC<-dplyr::select(trait.LIHC,"case_submitter_id","age_at_index","gender","ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stage","ajcc_pathologic_t")
trait.LIHC<-trait.LIHC[!duplicated(trait.LIHC),]
#
LIHC_count_used_id<-unlist(lapply(strsplit(colnames(TPM_log)[grep("TCGA",colnames(TPM_log))],split = "\\-"),function(x){paste0(x[1:3],collapse = "-")}))
summ_LIHC_count_used_id<-as.data.frame(table(LIHC_count_used_id))

trait.LIHC<-trait.LIHC[trait.LIHC$case_submitter_id%in%LIHC_count_used_id,]

for (i in trait.LIHC$case_submitter_id) {
  if(summ_LIHC_count_used_id[summ_LIHC_count_used_id$LIHC_count_used_id==i,2]>1){
    repno<-summ_LIHC_count_used_id[summ_LIHC_count_used_id$LIHC_count_used_id==i,"Freq"]
    trait.LIHC<-rbind(trait.LIHC,matrix(data = rep(unlist(trait.LIHC[trait.LIHC$case_submitter_id==i,]),repno-1),ncol = 7,nrow = repno-1,byrow = T,dimnames = list(NULL,colnames(trait.LIHC))))
  }
}
#sum(summ_LIHC_count_used_id$Freq)==nrow(trait.LIHC)
trait.LIHC<-trait.LIHC[match(LIHC_count_used_id,trait.LIHC$case_submitter_id),]
trait.LIHC$case_submitter_id<-colnames(TPM_log)[grep("TCGA",colnames(TPM_log))]
colnames(trait.LIHC)[1]<-"Samples"
#
phenotype<-left_join(x = phenotype,y = trait.LIHC,by="Samples")
phenotype<-phenotype[match(colnames(TPM_log),phenotype$Samples),]
phenotype<-unite(phenotype, "age", age,age_at_index,  remove = FALSE)
phenotype$age<-gsub(pattern = "NA_|_NA",replacement = "",phenotype$age)
phenotype$age<-as.numeric(phenotype$age)
phenotype<-unite(phenotype, "gender", gender.x,gender.y,  remove = FALSE)
phenotype$gender<-gsub(pattern = "NA_|_NA",replacement = "",phenotype$gender)
#
phenotype<-dplyr::select(phenotype,"Samples","Group","Batch","KI67","Villin","age","gender","ajcc_pathologic_m","ajcc_pathologic_n","ajcc_pathologic_stage","ajcc_pathologic_t","Vascular_tumor_thrombus")

#分期
phenotype$ajcc_pathologic_m<-mapvalues(x = phenotype$ajcc_pathologic_m,from = c("MX"),to = c(NA))
phenotype$ajcc_pathologic_n<-mapvalues(x = phenotype$ajcc_pathologic_n,from = c("NX","'--"),to = c(NA,NA))
phenotype$ajcc_pathologic_t<-mapvalues(x = phenotype$ajcc_pathologic_t,from = c("TX","'--"),to = c(NA,NA))
phenotype$ajcc_pathologic_stage<-mapvalues(x = phenotype$ajcc_pathologic_stage,from = c("'--"),to = c(NA))
#
rownames(phenotype)<-NULL
phenotype<-as.data.frame(phenotype)%>%column_to_rownames("Samples")
colnames(phenotype)[ncol(phenotype)]<-"VTT"

######five classfication
phenotype<-read.table("A://HCC/Bulkseq/WGCNA/HCC_mix/phenotype_onlyT.txt",sep = "\t")
meta<-read.table("A://HCC/Bulkseq/survival_class/HCC_mix_5.txt")
all.equal(rownames(phenotype),rownames(meta))
phenotype$catgory<-meta$catgory
write.table(x = phenotype,file = "A://HCC/Bulkseq/WGCNA/HCC_mix/phenotype_onlyT_5.txt",quote = F,sep = "\t",row.names = T,col.names = T)

