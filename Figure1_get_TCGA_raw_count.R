#TCGA RNA Seq 
library(rjson)
library(jsonlite)
library("rlist")
library("dplyr")
json_File <- fromJSON("A:/HCC/Bulkseq/TCGA_LIHC/TCGA_LIHC.json")
datasave<-list()
for( i in 1:length(json_File[["file_name"]])){
  datasave[[json_File[["file_name"]][i]]]<-json_File[["associated_entities"]][[i]]["entity_submitter_id"]
}
matchinfo<-list.rbind(datasave)
matchinfo<-cbind(rownames(matchinfo),matchinfo)
colnames(matchinfo)<-c("loc","TCGAID")
#
countsraw<-"A://HCC/Bulkseq/TCGA_LIHC/countsraw/"
countsloc<-paste0(countsraw,dir(countsraw,recursive = T))
counts_file<-lapply(matchinfo$TCGAID, function(x){
  filename<-matchinfo[which(matchinfo$TCGAID==x),"loc"]
  res<-data.table::fread(countsloc[grep(filename,countsloc)])
  res
})
names(counts_file)<-matchinfo$TCGAID
note<-"^unstranded$" # raw count 
# if the note is "tpm_unstranded" then get tpm
#过滤
counts_file_list<-lapply(names(counts_file), function(x){
  tmp<-counts_file[[x]][-(1:4),c(1,4,7)]
  tmp<-tmp[,c("gene_id","tpm_unstranded","unstranded")]
  colnames(tmp)[grep(note,colnames(tmp))]<-x
  tmp
})

names(counts_file_list)<-names(counts_file)

counts_merge<-matrix(data = NA,nrow = nrow(counts_file_list[[1]]),ncol = length(counts_file_list),dimnames = list(counts_file_list[[1]]$gene_id,names(counts_file)))

for (i in colnames(counts_merge)) {
  counts_merge[,i]<-as.numeric(unlist(as.data.frame(counts_file_list[[i]])[,grep("TCGA",colnames(counts_file_list[[i]]))]))
  print(i)
}
counts_merge<-as.data.frame(counts_merge)
#transfer
ids<-unlist(lapply(strsplit(rownames(counts_merge),"\\."),function(x){x[1]}))
counts_merge$"gene_ids"<-ids
#
load("A://HCC/Bulkseq/gtffile/Homo_sapiens.GRCh37.p13.gtf.RData")
gtf_used<-subset(GRCH37p13,type=="gene")
gtf_used<-dplyr::select(gtf_used,gene_id,gene_name)
gtf_used$gene_id<-unlist(lapply(strsplit(x = gtf_used$gene_id,split = "\\."),function(x){x[[1]]}))
#
counts_merge<-counts_merge[counts_merge$gene_ids%in%gtf_used$gene_id,]
counts_merge<-as.matrix(counts_merge)
rownames(counts_merge)<-gtf_used[match(counts_merge[,"gene_ids"],gtf_used$gene_id),"gene_name"]
counts_merge<-counts_merge[,-ncol(counts_merge)]
counts_merge_nu<-apply(counts_merge,2,as.numeric)
rownames(counts_merge_nu)<-rownames(counts_merge)
#QC
counts_merge_clean<-counts_merge_nu[rowSums(counts_merge_nu)>=10,]
index=order(rowMeans(counts_merge_clean),decreasing = T)
expr_ordered=counts_merge_clean[index,]
keep=!duplicated(rownames(expr_ordered))
expr_max=expr_ordered[keep,]
expr_max<-round(expr_max,digits = 2)
#
table(duplicated(rownames(expr_max)))
write.table(expr_max,"A://HCC/Bulkseq/TCGA_LIHC/expr_ma/rawcount/rawcount.txt",col.names = T,row.names = F,sep = "\t",quote = F)
expr_max <- data.table::fread("A://HCC/Bulkseq/TCGA_LIHC/expr_ma/rawcount/rawcount.txt")
writexl::write_xlsx(x = expr_max,path = "A://HCC/submit/0606submit/elife/sourcedataforfigure/Figure1_S1/data/TCGA_LIHC_raw_count.xlsx")
