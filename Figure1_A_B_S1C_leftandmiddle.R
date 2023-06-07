library(reshape2)
library("ggplot2")
library("rlist")
library("ggpubr")
#HCC_Mix
rawdata<-as.matrix(read.csv("A://HCC/Bulkseq/CIBERSORTx/output/singlecell/HCC_Mix_fraction.csv",header = T,row.names  = 1))
rawdata_ <- cbind(rownames(rawdata),rawdata)
colnames(rawdata_)[1]<-"Sample"
writexl::write_xlsx(x = as.data.frame(rawdata_),path = "A://HCC/submit/0606submit/elife/sourcedataforfigure/Figure1_S1/data/cell_fractions_raw.xlsx",col_names = T)
relativedata_raw<-apply(rawdata, 1,FUN = function(x){x[1:(ncol(rawdata)-4)]/x[ncol(rawdata)]})
relativedata<-reshape2::melt(relativedata_raw)
colnames(relativedata)<-c("celltype","sample","percent")
relativedata$celltype<-factor(relativedata$celltype,levels = rev(colnames(rawdata)[1:(ncol(rawdata)-4)]))
# for Mix
meta<-read.table("A://HCC/Bulkseq/CIBERSORTx/myownfiles/phen.txt",header = T)
only_T<-meta[meta$Group=="T","Samples"]
only_T_raw<-relativedata_raw[,colnames(relativedata_raw)%in%only_T]
only_T_forbar<-relativedata[relativedata$sample%in%only_T,]

only_T_forbar$"Absolute_TME_score"<-NA
for(i in unique(only_T_forbar$sample)){
  only_T_forbar[only_T_forbar$sample==i,"Absolute_TME_score"]<-only_T_forbar[only_T_forbar$sample==i,"percent"]*rawdata[i,ncol(rawdata)]
}

cell_sle<-list()
for (i in 1:ncol(only_T_raw)) {
  thre<-0.6
  tmp<-sort(only_T_raw[,i],decreasing = T)
  for (j in 1:length(tmp)) {
    if(sum(tmp[1:j])>thre){ 
      cell_sle[[colnames(only_T_raw)[i]]]<-names(tmp)[1:j]
      break}
  }
}
cell_sle_u<-names(sort(table(unlist(cell_sle)),decreasing = T)[1:5])

test<-only_T_raw[cell_sle_u,]
p<-pheatmap::pheatmap(test,show_colnames = T,scale = "column",cluster_rows = T,cutree_cols = 5)
mat_cluster <- test[p$tree_row$order, p$tree_col$order]
#Fig-S1A
only_T_forbar$sample<-factor(only_T_forbar$sample,levels = colnames(mat_cluster))
p1<-ggplot(only_T_forbar,aes(sample,Absolute_TME_score,fill=celltype))+geom_bar(stat="identity",position="stack")+theme_classic()+theme(axis.text.x = element_text(size=10,color = "black",face = "bold",angle = 45, hjust = 1, vjust = 1))
#ggsave(filename = "A://HCC/figs/bulk/bar/cellfraction_all.pdf",plot = p1,height = 10,width = 20)
#
test.clust <- rbind(test,cluster = cutree(p$tree_col, k = 5))
ann_col<-as.data.frame(test.clust[length(cell_sle_u)+1,])
colnames(ann_col)<-"catgory"
ann_col[ann_col$catgory==1,"catgory"]<-"pro_T"
ann_col[ann_col$catgory==2,"catgory"]<-"pro_meta"
ann_col[ann_col$catgory==3,"catgory"]<-"T_meta_mix"
ann_col[ann_col$catgory==4,"catgory"]<-"MemT_"
ann_col[ann_col$catgory==5,"catgory"]<-"NKC_"
ann_col$catgory<-factor(ann_col$catgory,levels = c("pro_T","T_meta_mix","pro_meta","NKC_","MemT_"))
labs.col <- colnames(test)
labs.col[grep("TCGA",labs.col)] <- ""
color_select<-function(class_){
  color<-as.data.frame(readxl::read_xlsx("A://HCC/配色.xlsx",sheet = "Sheet1"))[,1:2]
  colnames(color)<-c("clusters","color")  
  loc<-match(class_,color$clusters)
  target<-color[loc,"color"]
  names(target)<-color[loc,"clusters"]
  return(target)
}
p1<-pheatmap::pheatmap(test,show_colnames = T,scale = "column",cluster_rows = T,show_rownames = T
                       ,cutree_cols = 5,annotation_col = ann_col,labels_col = labs.col,cellwidth = 2.5,cellheight = 25,border_color = "grey",fontsize = 12,annotation_colors = list(catgory=color_select(class_ = levels(ann_col$catgory))))
#write.table(ann_col,"A://HCC/Bulkseq/survival_class/HCC_mix_5.txt",sep = "\t",quote = F,row.names = T,col.names = T)


#####t-test
library(ggplot2)
library(tidyverse)
for_t_test<-only_T_forbar
meta<-read.table("A://HCC/Bulkseq/survival_class/HCC_mix_5.txt")
for_t_test$"TB"<-mapvalues(x = for_t_test$sample,from = rownames(meta),to = meta$catgory)
P_result<-list()
for(i in levels(for_t_test$celltype)){
  cellper<-subset(for_t_test,celltype==i)
  tmp<-t.test(cellper[cellper$TB=="pro_T","percent"],alternative = "two.sided",cellper[cellper$TB=="pro_meta","percent"])
  #
  P_result[[i]][["p.value"]]<-tmp$p.value
  P_result[[i]][["estimate"]]<-tmp$estimate
}
P_result_p<-list.rbind(lapply(P_result,function(x){x$p.value}))
P_result_p_adj<-p.adjust (P_result_p, method="fdr")
P_result_p_adj<-cbind(rownames(P_result_p),P_result_p_adj)
colnames(P_result_p_adj)<-c("cell","p.adj")
P_result_sig<-as.data.frame(na.omit(P_result_p_adj[as.numeric(P_result_p_adj[,2])<=0.05,]))

allpic<-list()
getbox<-function(cell){
  forbox<-for_t_test[for_t_test$celltype==cell&!for_t_test$TB%in%c("NKC_","MemT_"),]
  forbox$celltype<-as.character(forbox$celltype)
  forbox$status<-factor(forbox$TB,levels = c("pro_T","T_meta_mix","pro_meta","NKC_","MemT_"))
  forbox$percent<-forbox$percent*100
  #,label = c("p.signif")
  p<-ggboxplot(forbox, x="status", y="percent", color = "status", title = cell,palette = color_select(class_ = levels(forbox$status)),font.label = list(size = 15, color = "black"))#+scale_y_continuous(limits=c(0,ceiling(max(forbox$percent))),breaks=round(seq(0,ceiling(max(forbox$percent)),length.out=5),2))
  #my_comparisons <- list(c("N","T"))
  #my_comparisons <- list(c("TBN_T_noJKC","TBY_T_noJKC"))
  my_comparisons <- list(c("pro_T","T_meta_mix"),c("T_meta_mix","pro_meta"),c("pro_meta","pro_T"))
  p1<-p+stat_compare_means(comparisons = my_comparisons, method = "t.test",size=8,paired = F,method.args = list(alternative="two.sided"),vjust = 0.6)+theme(text = element_text(size = 8))
  allpic[[cell]]<<-p1
  #ggsave(paste0("A://HCC/figs/bulk/box/cellfraction_TCGA/",cell,".pdf"),p1,height = 7,width = 8)
}
for(i in P_result_sig$cell){
  getbox(cell = i)}

up<-c()
dn<-c()
lapply(names(P_result),function(x){
  ifelse(P_result[[x]][["estimate"]][1]<P_result[[x]][["estimate"]][2],up<<-c(up,x),dn<<-c(dn,x))})


disrow<-2
p_up<-ggarrange(plotlist = allpic[names(allpic)%in%up],nrow = disrow,ncol = ceiling(length(allpic[names(allpic)%in%up])/disrow))
p_dn<-ggarrange(plotlist = allpic[names(allpic)%in%dn],nrow = disrow,ncol = ceiling(length(allpic[names(allpic)%in%dn])/disrow))

### a wide-type plot, The significance  symbols are manually added based on p_up and p_dn

cells<-read.table("A://HCC/Bulkseq/survival_class/cells.txt",header = T)
library("plyr")
forbox_wide<-for_t_test
colnames(forbox_wide)[3]<-"Percentage(%)"
forbox_wide$`Percentage(%)`<-100*forbox_wide$`Percentage(%)`
forbox_wide<-forbox_wide[(!forbox_wide$celltype%in%c("patient.specific.macropahge","proliferative.macrophage","proliferative.CD4.T","proliferative.CD8.T"))&forbox_wide$celltype%in%cells$cells,]
forbox_wide<-forbox_wide[!forbox_wide$TB%in%c("MemT_","NKC_"),]
forbox_wide$TB<-factor(forbox_wide$TB,levels = c("pro_T","T_meta_mix","pro_meta"))
forbox_wide$celltype<-factor(forbox_wide$celltype,levels = cells$cells)
p1<-ggplot(data = forbox_wide,aes(x=celltype,y = `Percentage(%)`,fill=TB))+stat_boxplot(geom = "errorbar",width=0.3,position = position_dodge(0.8))+geom_boxplot(width=0.5,position = position_dodge(0.8))+theme_classic()+theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))+scale_fill_manual(values = color_select(class_ = levels(forbox_wide$TB)))

## Fig S1B left and middle
#
ann_col<-cbind(rownames(ann_col),ann_col)
colnames(ann_col)[1]<-"SampleID"

top5cells<-test
allcells<-only_T_raw
#
forpca_data<-allcells # or forpca_data<-top5cells

pca<-prcomp(t(forpca_data),scale.=F,center=F)
library("tibble")
library("dplyr")
score <- inner_join(pca$x %>% data.frame() %>%
                      rownames_to_column("SampleID") %>%
                      select(c(1:3)),
                    ann_col,
                    by="SampleID")
p1<-ggplot(score,aes(x=PC1,y=PC2))+
  geom_point(aes(color=score[,"catgory"]),size=1)+scale_color_manual(values=color_select(class_ = score$catgory),name="catgory")+stat_ellipse(aes(fill=score[,"catgory"]),type="norm",geom="polygon",alpha=0.2,color=NA)+scale_fill_manual(values=color_select(class_ = levels(score$catgory)))+theme_classic()
#ggsave(filename = "A://HCC/figs/bulk/PCA/allcells.pdf",p1,width = 7,height = 6)

