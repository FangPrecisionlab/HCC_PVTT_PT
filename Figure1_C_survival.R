#----###################fromTCGA-----###
library(tinyarray)
library("survival")
library("survminer")
expr_max<-read.table("A://HCC/Bulkseq/VIPER/preforARACNe/TPM/log/TCGA_LIHC_allclean_log2.txt",header = T,row.names = 1,check.names = F)
expr_max_<-cbind(rownames(expr_max),expr_max)
colnames(expr_max_)[1]<-"gene"
writexl::write_xlsx(x = expr_max_,path = "A://HCC/submit/0606submit/elife/sourcedataforfigure/Figure1_S1/data/TCGA_LIHC_log2TPM.xlsx")
#
Normal<-colnames(expr_max)[unlist(lapply(strsplit(colnames(expr_max),split = "-"),function(x){substr(x[4],1,2)}%in%c(10:29)))]
Tumor<-colnames(expr_max)[!unlist(lapply(strsplit(colnames(expr_max),split = "-"),function(x){substr(x[4],1,2)}%in%c(10:29)))]
expr_T <- expr_max[,make_tcga_group(expr_max)=="tumor"]
#
clinical<-read.table("A://HCC/Bulkseq/TCGA_LIHC/clinical/clinical.tsv",header = T,check.names = F,sep = "\t",quote = "")
t_needed=c("case_submitter_id",
           "vital_status",
           "days_to_last_follow_up",
           "days_to_death",
           "age_at_index",
           "ajcc_pathologic_t")
clinical<-clinical[,t_needed]
clinical=clinical[clinical$vital_status %in% c('Alive','Dead'),]
#
clinical$days_to_last_follow_up[clinical$days_to_last_follow_up=="'--"] = 0 #is.na()用于返回是否为缺失值
clinical$days_to_death[clinical$days_to_death=="'--"] = 0   
clinical$days<-ifelse(clinical$vital_status=='Alive',clinical$days_to_last_follow_up,clinical$days_to_death)
clinical$days<-as.numeric(clinical$days)
clinical$month=round(clinical$days/30,2) #以month为单位，保留两位小数
clinical<-clinical[!clinical$age_at_index=="'--",]
clinical$age_at_index<-as.numeric(clinical$age_at_index)
clinical$age_group = ifelse(clinical$age_at_index>median(clinical$age_at_index),'old','young')
clinical$"sample"<-"blank"
clinical<-clinical[!duplicated(clinical$case_submitter_id),]
#write.table(clinical,"A://HCC/Bulkseq/TCGA_LIHC/clinical/clinical_filtered.tsv")
#
clinical<-read.table("A://HCC/Bulkseq/TCGA_LIHC/clinical/clinical_filtered.tsv")
meta<-clinical
matched<-colnames(expr_T)[unlist(lapply(unique(meta$case_submitter_id), function(x){grep(x,colnames(expr_T))}))]

allpatient<-unlist(lapply(strsplit(matched,split = "\\-"),FUN = function(x){paste0(x[1:3],collapse = "-")}))

multi<-allpatient[duplicated(allpatient)]
single<-allpatient[!allpatient%in%multi]
meta<-meta[meta$case_submitter_id%in%single,]

meta[,"sample"]<-colnames(expr_T)[unlist(lapply(single, function(x){grep(x,colnames(expr_T))}))]
#表达矩阵中的重复样本
ins<-colnames(expr_T)[unlist(lapply(unique(multi), function(x){grep(x,colnames(expr_T))}))]
#
app<-matrix(data = NA,nrow = length(ins),ncol = ncol(meta),dimnames = list(NULL,colnames(meta)))
app[,"sample"]<-ins
app[,"case_submitter_id"]<-unlist(lapply(strsplit(app[,"sample"],split = "\\-"),FUN = function(x){paste0(x[1:3],collapse = "-")}))

library("rlist")
app<-list.rbind(apply(app, MARGIN = 1,function(x){
  x[2:9]<-clinical[grep(x[1],clinical$case_submitter_id)[1],colnames(app)][2:9]
  x
}))
#
meta2<-as.data.frame(lapply(rbind(meta,app),unlist))
metaforcox<-dplyr::select(meta2,"sample","vital_status","case_submitter_id","month")
colnames(metaforcox)<-c("sample","event","X_PATIENT","time")
metaforcox[metaforcox$event=="Alive","event"]<-0
metaforcox[metaforcox$event=="Dead","event"]<-1
metaforcox<-as.data.frame(lapply(metaforcox, unlist))
#write.table(metaforcox,"A://HCC/Bulkseq/TCGA_LIHC/clinical/metaforcox.tsv")

###
metaforcox<-read.table("A://HCC/Bulkseq/TCGA_LIHC/clinical/metaforcox.tsv")
###
expr_T_dn<-expr_T[,metaforcox$sample]
expr_T_dn = expr_T_dn[apply(expr_T_dn, 1, function(x){sum(x==0) < 0.25 *ncol(expr_T_dn)}),]
#five classes groups
classes<-read.table("A://HCC/Bulkseq/survival_class/HCC_mix_5.txt")
classes<-cbind(row.names(classes),classes)
colnames(classes)<-c("ID","group2")
classes<-classes[metaforcox$sample,]
table(classes$ID==metaforcox$sample)
metaforcox$group2<-classes$group2
metaforcox$event<-as.integer(metaforcox$event)
#
survData = Surv(time = metaforcox$time,                 #生存时间数据
                event = metaforcox$event) #判断结局，完全数据/截尾数据

KMfit <- survfit(survData ~ metaforcox$group2)  # ~ 后是指定的分组
##画图
color_select<-function(class_){
  color<-as.data.frame(readxl::read_xlsx("A://HCC/配色.xlsx",sheet = "Sheet1"))[,1:2]
  colnames(color)<-c("clusters","color")  
  loc<-match(class_,color$clusters)
  target<-color[loc,"color"]
  names(target)<-color[loc,"clusters"]
  return(target)
}
metaforcox$group2<-factor(metaforcox$group2,levels = c("pro_T","T_meta_mix","pro_meta","NKC_","MemT_"))
#pdf(file = "A://HCC/figs/bulk/survial/survival_5classes.pdf",width = 12,height = 10.5)
ggsurvplot(KMfit,                     #拟合对象
           data = metaforcox,               #变量数据来源
           pval = TRUE,               #P值
           surv.median.line = "hv",   #中位生存时间线
           risk.table = TRUE,         #风险表
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
           palette = unname(color_select(class_ = levels(metaforcox$group2))),
           ggtheme = theme_classic()
)
dev.off()
pairwise_survdiff(Surv(time, event) ~ group2,
                  data = metaforcox)
