library("viper")
library("aracne.networks")
library("org.Hs.eg.db")
library(AnnotationDbi)
library("minet")
library("rlist")
library("data.table")
library(GenomicDataCommons)

###start viper###
###import the files
#
network_file<-data.table::fread("A://HCC/Bulkseq/VIPER/downstream/TCGA/tf/consoli/network.txt",header = T)
adj_file<-as.matrix(network_file[,-"pvalue"])
write.table(adj_file,"A://HCC/Bulkseq/VIPER/downstream/TCGA/tf/consoli/network.adj",col.names = F,row.names = F,quote = F,sep = "\t")
##mak regulon obj
regul <- aracne2regulon(afile = "A://HCC/Bulkseq/VIPER/downstream/TCGA/tf/consoli/network.adj", exp, verbose = T,format = c("3col"))
save(list = "regul",file = "A://HCC/Bulkseq/VIPER/downstream/TCGA/tf/consoli/TCGA_regul_tf.rda")


############start
exp<-as.matrix(read.table("A://HCC/Bulkseq/VIPER/preforARACNe/TPM/TCGA_LIHC_allclean.txt",header = T,check.names = F,row.names = 1))
load("A://HCC/Bulkseq/VIPER/downstream/TCGA/tf/consoli/TCGA_regul_tf.rda")
print(regul)
###group
#T_N
Normal<-colnames(exp)[unlist(lapply(strsplit(colnames(exp),split = "-"),function(x){substr(x[4],1,2)}%in%c(10:29)))]
Tumor<-colnames(exp)[!unlist(lapply(strsplit(colnames(exp),split = "-"),function(x){substr(x[4],1,2)}%in%c(10:29)))]
###5分类法
classes<-read.table("A://HCC/Bulkseq/survival_class/HCC_mix_5.txt")
colnames(classes)<-c("ID","group2")
#
pro_T<-classes[classes$group2=="pro_T","ID"]
T_meta_mix<-classes[classes$group2=="T_meta_mix","ID"]
pro_meta<-classes[classes$group2=="pro_meta","ID"]
MemT_<-classes[classes$group2=="MemT_","ID"]
NKC_<-classes[classes$group2=="NKC_","ID"]

###mak signature
signature <- rowTtest(x = exp[,colnames(exp)%in%pro_meta], y = exp[,colnames(exp)%in%pro_T])
signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) *sign(signature$statistic))[, 1]
#
nullmodel<-ttestNull(x = exp[,colnames(exp)%in%pro_meta], y = exp[,colnames(exp)%in%pro_T],per = 1000,repos = TRUE, verbose = T)
save("nullmodel",file = "A://HCC/Bulkseq/VIPER/requireData/5classes/nullmodel.RData")
###msVIPER
mrs <- msviper(signature, regul, nullmodel, verbose = T)
summary(mrs)
pdf(file = "A://HCC/figs/bulk/viper/5classes/noboot_LIHC_pro_meta_T.pdf",width = 10,height = 10)
plot(mrs, mrs=12,cex = .9)
dev.off()
######Leading-edge analysis
mrs <- ledge(mrs)
save("signature","mrs",file = "A://HCC/Bulkseq/VIPER/requireData/5classes/nobootmrs_metavsT.RData")

#####Bootstrap msVIPER,采用重采样技术可以减少基因表达特征上的异常样本

signature <- bootstrapTtest(x=exp[,colnames(exp)%in%pro_meta], y = exp[,colnames(exp)%in%pro_T], verbose = T)
mrs <- msviper(signature, regul, nullmodel, verbose = T)
mrs <- bootstrapmsviper(mrs, "mode")
save("signature","mrs",file = "A://HCC/Bulkseq/VIPER/requireData/5classes/bootmrs.RData")
load("A://HCC/Bulkseq/VIPER/requireData/5classes/bootmrs.RData")
pdf(file = "A://HCC/figs/bulk/viper/5classes/boot_LIHC_pro_meta_T.pdf",width = 10,height = 10)
plot(mrs, cex = .9)
dev.off()
pdf(file = "A://HCC/figs/bulk/viper/5classes/boot_LIHC_pro_meta_T_TOP30.pdf",width = 10,height = 10)
plot(mrs, 30,cex = .9)
dev.off()
