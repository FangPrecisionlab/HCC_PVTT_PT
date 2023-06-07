library(maftools)
library(mclust)
library(NMF)
library(pheatmap)
library(barplot3d)
library(openxlsx)
#tcga
file_LIHC<-"A://HCC/Wesseq/TCGAWXS/all_rebarcode.maf"
tcga_clin<-"A://HCC/Wesseq/TCGAWXS/clinical_maf.tsv"
tcga_maf<-read.maf(maf = file_LIHC,clinicalData = tcga_clin)
#local
file_local<-"A://HCC/Wesseq/finalresult/maf/all.maf"
local_maf<-read.maf(file_local)
LIHC_MIX_maf<-merge_mafs(mafs = list(tcga_maf,local_maf))
###分组
meta<-read.table("A://HCC/Bulkseq/survival_class/HCC_mix_5.txt")
sampleid<-list()
for (i in unique(meta$catgory)) {
  local_samples<-meta$SampleID[grep("^(?!TCGA)",meta$SampleID,perl = T)]
  tmp<-rownames(meta)[which(meta$catgory==i)]
  #
  local<-intersect(tmp,local_samples)
  TCGA<-tmp[!tmp%in%local]
  #
  sampleid[[i]]<-c(local,unlist(lapply(strsplit(TCGA,"-"),function(x){paste0(x[1:3],collapse = "-")})))
}
#
sample_maf<-list()
for (i in unique(meta$catgory)) {
  sample_maf[[i]]<- subsetMaf(maf = LIHC_MIX_maf,tsb = sampleid[[i]])
}
Metas.vs.pro_T <- mafCompare( m1 = sample_maf$pro_meta, m2 = sample_maf$pro_T,m1Name = 'Metas', m2Name = 'pro_T', minMut = 5)
Metas.vs.Mix<-mafCompare( m1 = sample_maf$pro_meta, m2 = sample_maf$T_meta_mix,m1Name = 'Metas', m2Name = 'Mix_maf', minMut = 5)
pro_T.vs.Mix<-mafCompare( m1 = sample_maf$pro_T, m2 = sample_maf$T_meta_mix,m1Name = 'pro_T', m2Name = 'Mix_maf', minMut = 5)
#
forestPlot(mafCompareRes = Metas.vs.pro_T, pVal = 0.01, color = c('#1DAA90','#FF4500'), geneFontSize = 1,lineWidth = 2,titleSize = 2)
forestPlot(mafCompareRes = Metas.vs.Mix, pVal = 0.05, color = c('#FFB6C1','#FF4500'), geneFontSize = 1,lineWidth = 2,titleSize = 2)
forestPlot(mafCompareRes = pro_T.vs.Mix, pVal = 0.05, color = c('#FFB6C1','#1DAA90'), geneFontSize = 1,lineWidth = 2,titleSize = 2)
#
M1.pfam = pfamDomains(maf = sample_maf$pro_meta, top = 10)
M2.pfam = pfamDomains(maf = sample_maf$pro_T, top = 10)
