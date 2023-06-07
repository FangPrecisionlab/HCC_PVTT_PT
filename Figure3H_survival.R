library(tinyarray)
library("survival")
library("survminer")
library("data.table")
library("tibble")
library("rlist")
#
expr<-fread("A://HCC/Bulkseq/VIPER/preforARACNe/TPM/log/TCGA_LIHC_allclean_log2.txt",check.names = F)
expr<-column_to_rownames(.data = expr,var = "gene")
#
Normal<-colnames(expr)[unlist(lapply(strsplit(colnames(expr),split = "-"),function(x){substr(x[4],1,2)}%in%c(10:29)))]
Tumor<-colnames(expr)[!unlist(lapply(strsplit(colnames(expr),split = "-"),function(x){substr(x[4],1,2)}%in%c(10:29)))]
#
expr_T <- expr[,make_tcga_group(expr)=="tumor"]
metaforcox<-read.table("A://HCC/Bulkseq/TCGA_LIHC/clinical/metaforcox.tsv")
###
expr_T<-expr_T[,metaforcox$sample]
expr_T = expr_T[apply(expr_T, 1, function(x){sum(x==0) < 0.25 *ncol(expr_T)}),]
#
survData = Surv(time = metaforcox$time,                 #生存时间数据
                event = metaforcox$event) #判断结局，完全数据/截尾数据

#分组,按照基因表达来分

getsurv<-function(genes){
  expr_T_dn = expr_T[rownames(expr_T)%in%genes,]
  rank_res<-as.data.frame(sort(apply(as.matrix(expr_T_dn), 2, mean),decreasing = T))
  colnames(rank_res)<-"mean_value"
  #
  line<-fivenum(rank_res$mean_value)
  high<-rownames(rank_res)[rank_res$mean_value>=line[3]]
  low<-rownames(rank_res)[rank_res$mean_value<line[3]]
  #
  metaforcox$group2<-NA
  metaforcox[metaforcox$sample%in%high,"group2"]<-"high"
  metaforcox[metaforcox$sample%in%low,"group2"]<-"low"
  
  
  KMfit <- survfit(survData ~ metaforcox$group2)  # ~ 后是指定的分组
  p_value<<-pairwise_survdiff(Surv(time, event) ~ group2,data = metaforcox)$p.value
  Exp_high_mid<-surv_median(KMfit)[grep("high",surv_median(KMfit)$strata),"median"]
  Exp_low_mid<-surv_median(KMfit)[grep("low",surv_median(KMfit)$strata),"median"]
  Group_info_all[Group_info_all$Gene==genes,c("Pvalue","Exp_high_mid","Exp_low_mid")]<<-c(p_value,Exp_high_mid,Exp_low_mid)
  
  
  p<-ggsurvplot(KMfit,                     #拟合对象
                data = metaforcox,               #变量数据来源
                pval = F,               #P值
                surv.median.line = "hv",   #中位生存时间线
                risk.table = F,         #风险表
                xlab = "Follow up time(m)",  #x轴标签
                break.x.by = 10,          #x轴刻度间距
                #legend = c(0.8,0.75),     #图例位置
                #legend.title = "",        #图例标题
                #legend.labs = c("old", "young"),  #图例分组标签
                font.x = c(12, "bold", "red"),
                font.y = c(12, "bold", "red"),
                font.tickslab = c(12, "bold", "black"),
                font.legend =  c(12, "bold", "black"),
                censor.size=7,
                # palette = unname(color_select(class_ = levels(metaforcox$group2))),
                ggtheme = theme_classic()
  )
  p2<-p$plot+annotate(geom="text", x=fivenum(metaforcox$time)[2], y=0.2, label=paste0(genes,"\n","p = ",round(p_value,4)),colour="red", size=5)+annotate(geom="text", x=surv_median(KMfit)$median[1], y=0, label=surv_median(KMfit)$median[1],colour="#F08080", size=5)+annotate(geom="text", x=surv_median(KMfit)$median[2], y=0, label=surv_median(KMfit)$median[2],colour="#008080", size=5)
  
  #ggsave(filename = paste0("A://HCC/figs/bulk/survial/genes/PPI_WGCNA/",genes,".pdf"),plot = p2,width = 7,height = 6.5)
}
##
library("dplyr")
library("plyr")
geneslist_dir<-dir("A:/HCC/Bulkseq/STRING",full.names = T,recursive = T)
geneslist_dir<-geneslist_dir[grep("anncol",geneslist_dir)]
geneslist<-list()
for (i in geneslist_dir) {
  info<-unlist(strsplit(x=i,split = "\\/",perl = T))
  geneslist[[info[5]]][[unlist(strsplit(info[6],"_"))[1]]]<-read.csv(i)
}

geneslist<-lapply(X = geneslist,FUN = function(x){lapply(x,function(y){y[!y$font_color=="#000000",]})})

Group_info_all<-list()
lapply(names(geneslist), function(x){
  Group_info<-list.rbind(lapply(geneslist[[x]],function(x){x}))
  Group_info<-cbind(rownames(Group_info),x,Group_info)
  colnames(Group_info)[c(1,2,3,5)]<-c("Cluster","Group","Gene","Module")
  Group_info<-select(Group_info,c("Cluster","Group","Gene","Module"))
  Group_info$Module<-mapvalues(x = Group_info$Module,from = c("#800000","#556B2F","#4B0082","#00CED1"),to = c("black","black","blue","blue"))
  Group_info_all[[x]]<<-Group_info
}
)
Group_info_all<-as.data.frame(list.rbind(Group_info_all))
Group_info_all[,c("Pvalue","Exp_high_mid","Exp_low_mid")]<-NA

for (i in Group_info_all$Gene) {
  getsurv(i)
}

Group_info_all$Cluster<-unlist(lapply(strsplit(Group_info_all$Cluster,split = "\\."),function(x){x[1]}))

#write.table(Group_info_all,"A://HCC/Bulkseq/survial/WGCNA_PPI.txt",quote = F,sep = "\t",col.names = T,row.names = F)

Metas_sur<-Group_info_all[Group_info_all$Group=="Metas",]
table(Metas_sur$Pvalue<0.05&Metas_sur$Exp_high_mid<Metas_sur$Exp_low_mid)


Pro_sur<-Group_info_all[Group_info_all$Group=="ProT",]
table(Pro_sur$Pvalue<0.05&Pro_sur$Exp_high_mid>Pro_sur$Exp_low_mid)
