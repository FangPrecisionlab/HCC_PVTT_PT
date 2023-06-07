#加载库
library(MASS)
library(class)
library(cluster)
library(impute)
library(Hmisc)
library(WGCNA)
library("DESeq2")
library("data.table")
library("tibble")
###############第零部分，读取数据#############
#参数设置
options(stringsAsFactors = F)# 在读入数据时，遇到字符串后，将其转换成因子，连续型变量要改为FALSE #工作目录
#setwd(dir = "A://HCC/Bulkseq/WGCNA/MADTOP5000/")
setwd(dir = "A://HCC/Bulkseq/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/")
#读取数据，用tpm,针对样本用hclust,所以转置一下,行为样本，列为基因
#TPM_log<-read.table("A://HCC/Bulkseq/degresult/TPM_log2.txt",header = T,row.names = 1,check.names = F)
TPM_log<-fread("A://HCC/Bulkseq/WGCNA/HCC_mix/mix_TPM_T_log2.txt",header = T,check.names = F)
TPM_log<-TPM_log%>%column_to_rownames("gene")
datExpr_raw<-as.data.frame(TPM_log,check.names=F)
dim(datExpr_raw)
###############第一部分，数据筛选#############
#筛选中位数绝对偏差（MAD）前75%的基因，至少MAD要大于0.01，并且不能用差异表达基因
#m.mad <- apply(datExpr_raw,1,mad)
#方式1：中位数偏差前75%的基因,并且还得大于0.01
#datExpr <- datExpr_raw[which(m.mad>max(quantile(m.mad,probs=seq(0,1,0.25))[2],0.01)),]
#方式2：MAD大于0.01
#datExpr <- datExpr_raw[which(m.mad>0.01),]
#方式3：mad前5000个基因
datExpr <- datExpr_raw[order(apply(datExpr_raw,1,mad), decreasing = T)[1:5000],]
#转换为样品在行，基因在列的矩阵
datExpr <- t(datExpr)
#check 检查缺失值
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
#QC 如果gsg$allOK的结果为TRUE，证明没有缺失值，可以直接下一步
#如果gsg$allOK的结果为TRUE，证明没有缺失值，可以直接下一步
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(datExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(datExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr)
#识别离群值（异常值）
sampleTree = hclust(dist(datExpr), method = "average")
pdf("A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/sample_clus.pdf",height = 20,width = 70)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
#去除离群样本
abline(h = 170, col = "red")  #划定需要剪切的枝长
clust = cutreeStatic(sampleTree, cutHeight = 170, minSize = 10)
table(clust)
keepSamples = (clust==1)  #保留非离群(clust==1)的样本
datExpr = datExpr[keepSamples, ]  #去除离群值后的数据，TOP5000不用截断
#再画一下
sampleTree = hclust(dist(datExpr), method = "average")
pdf("A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/sample_clus_rem.pdf",height = 20,width = 70)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
#表型数据
#traitData = xlsx::read.xlsx("A://HCC/Seqinfo/clinic_WGCNA.xlsx",sheetIndex = 1,header = T)
#traitData = read.table("A://HCC/Bulkseq/WGCNA/HCC_mix/phenotype_onlyT.txt",sep = "\t")
traitData = read.table("A://HCC/Bulkseq/WGCNA/HCC_mix/phenotype_onlyT_5.txt",sep = "\t")

#for local
#rownames(traitData) = traitData[,1]
#traitData<-traitData[, -1]

#for mixT
traitData<-traitData[rownames(datExpr),]
traitData<-traitData[,!unlist(lapply(apply(traitData, 2, unique),function(x){length(x)==1}))]
#用颜色代表关联度
collectGarbage()
#可视化表型数据与基因表达量数据的联系，重构样本聚类树
sampleTree2 = hclust(dist(datExpr), method = "average")
#
disper<-c("Batch","gender","ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m","ajcc_pathologic_stage","VTT","catgory")
#
traitData$Batch<-factor(x = traitData$Batch,levels = c("LOCAL","LIHC"))
traitData$gender<-factor(x = traitData$gender,levels = c("female","male"))
traitData$ajcc_pathologic_t<-factor(x = traitData$ajcc_pathologic_t,levels = c("T1","T2","T2a","T2b","T3","T3a","T3b","T4"))
traitData$ajcc_pathologic_n<-factor(x = traitData$ajcc_pathologic_n,levels = c("N0","N1"))
traitData$ajcc_pathologic_m<-factor(x = traitData$ajcc_pathologic_m,levels = c("M0","M1"))
traitData$ajcc_pathologic_stage<-factor(x = traitData$ajcc_pathologic_stage,levels = c("Stage I","Stage II","Stage III","Stage IIIA","Stage IIIB","Stage IIIC","Stage IV","Stage IVA","Stage IVB"))
traitData$VTT<-factor(x = traitData$VTT,levels = c(0,1))
traitData$catgory<-factor(x = traitData$catgory,levels = c("pro_T","T_meta_mix","pro_meta","NKC_","MemT_"))
#
traitData_nu<-traitData
for(i in disper){
  traitData_nu[,i]<-as.numeric(factor(traitData_nu[,i]))
}
#
traitColors = numbers2colors(traitData_nu, signed = T) #用颜色代表关联度
set.seed(1)
pdf(file = "A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/sample_clust_trait.pdf",width = 70,height = 40)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitData_nu),
                    main = "Sample dendrogram and trait heatmap",marAll = c(4,10,3,4))
dev.off()
#图片结果解释了临床数据和基因表达量的关联程度，保存数据
getwd()
save(datExpr, traitData, traitData_nu,file = "HCC_local_dataInput.RData")

#############第二部分筛选软阈值###########
type = "unsigned"
powers1 <- c(seq(1, 10, by=1), seq(12, 28, by=2))
sft = pickSoftThreshold(datExpr, powerVector=powers1, 
                        networkType=type, verbose=5,RsquaredCut = 0.9)
#绘图筛选
pdf("A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/sft.pdf",width = 15,height = 10)
par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers1,cex=cex1,col="red")
# 筛选标准。R-square=0.9 第一个高于0.9的
abline(h=0.90,col="red")
# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers1, 
     cex=cex1, col="red")
dev.off()
power = sft$powerEstimate
#经验power,适用于样本数少的情况
#也有说50左右的连通度，20左右的中位数
nSamples<-nrow(datExpr)
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))       
                 )
  )
}
power
rm(list = ls())

################################第三部分，构建加权共表达网络###################################
##在Rstudio中重启R，然后library(WGCNA)
#search()
##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值,9
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
# goto local linux
setwd("/mnt/hgfs/Win_Cent/WGCNA_HCC_mix_onlyT_h2/")
load("HCC_local_dataInput.RData")
#设置分析参数
power=8
MCH<-0.2
type = "unsigned"
corType = "pearson"
power = power
nGenes = ncol(datExpr)
objname <- "HCC_mix_T_h2"
#开启多线程
enableWGCNAThreads()
#
net = blockwiseModules(datExpr, #
                       power = power, #
                       maxBlockSize = nGenes,#
                       TOMType = type, #
                       minModuleSize = 30,#
                       reassignThreshold = 0, 
                       mergeCutHeight = MCH,#
                       numericLabels = TRUE, #
                       pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, #
                       corType = corType, #
                       networkType = type,#
                       loadTOMs=TRUE,
                       saveTOMFileBase = paste0(objname, ".tom"),
                       verbose = 3,#
                       nThreads=11)
save(net,file="./net.RData")

###now back to locl PC 
load("./net.RData")
#查看划分的模块数和每个模块里面包含的基因个数
table(net$colors)#模块0是无法识别的基因数,对应灰色块，未分配到任何模块的基因
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
moduleColors = labels2colors(net$colors)
sort(table(moduleColors),decreasing = T)
moduleLabels = net$colors
#如果需要修改树、模块成员、和模块合并标准，该包的recutBlockwiseTrees函数可以对结果进行修改，而无需重复计算网络和树状图。（推荐用第二种方法分步法实现）
pdf("A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/Dendrogram.pdf",height = 10,width = 10)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()



##################第四部分，关联基因模块与表型##################
load("./HCC_local_dataInput.RData")
load("./net.RData")
#####模块-表型数据关联并识别重要基因
#注意：以下代码专门针对转移和致瘤表型
getheat<-function(disp_var){
library(stringr)
traitData[,disp_var] <- as.factor(traitData[,disp_var])
options(na.action='na.pass')
design <- model.matrix(~0+traitData[,disp_var])
if(nrow(design))
colnames(design) <- levels(traitData[,disp_var]) #get the group
#
#method1
MES0<-net$MEs
colnames(MES0) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MES0),"ME",""))))
#method2 小数位数影响order结果
# moduleColors <- labels2colors(net$colors)
# MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  #Calculate module eigengenes.
MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other
modTraitCor <- cor(MEs,design,use = "p")
nSamples<-ncol(datExpr)
modTraitP <- corPvalueStudent(modTraitCor,nSamples)
textMatrix <- paste0(signif(modTraitCor,2),"\n(",
                     signif(modTraitP,1),")")
dim(textMatrix) <- dim(modTraitP)
pdf("A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/M_T_rela_cat.pdf",height = 8,width = 8)
par(mar=c(5, 9, 3, 3))
p1<-labeledHeatmap(Matrix = modTraitCor, 
               xLabels = colnames(design), 
               yLabels = colnames(MEs), 
               cex.lab = 1, 
               ySymbols = colnames(MEs), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, 
               setStdMargins = F, 
               cex.text = 0.7, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
return(design)
}
design<-getheat(disp_var = "catgory")

### 模块与表型的相关性boxplot图 
library(gplots)
library(ggpubr)
library(grid)
library(gridExtra) 
library("stringr")
#
color_select<-function(class_){
  color<-as.data.frame(readxl::read_xlsx("A://HCC/配色.xlsx",sheet = "Sheet1"))[,1:2]
  colnames(color)<-c("clusters","color")  
  loc<-match(class_,color$clusters)
  target<-color[loc,"color"]
  names(target)<-color[loc,"clusters"]
  return(target)
}
MES0<-net$MEs
colnames(MES0) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MES0),"ME",""))))
MEs <- orderMEs(MES0) 
mes_group <- merge(MEs,traitData,by="row.names") 
mes_group<-subset(mes_group,!catgory%in%c("NKC_","MemT_"))
mes_group$catgory<-factor(mes_group$catgory,levels = c("pro_T","T_meta_mix","pro_meta"))
draw_ggboxplot <- function(data,Module="Module",group="group"){
  ggboxplot(data,x=group, y=Module,
            ylab = paste0(Module),
            xlab = group,
            fill = group,
            #palette = "jco",
            #add="jitter",
            legend = "")+scale_fill_manual(values=color_select(class_ = levels(data$catgory)),name="catgory") +stat_compare_means(comparisons = list(c("pro_T","T_meta_mix"),c("T_meta_mix","pro_meta"),c("pro_T","pro_meta")),method = "t.test")}
p<-list()
for(i in unique(labels2colors(net$colors))){
 p[[i]]<- draw_ggboxplot(mes_group,  group = "catgory",Module = paste0("ME",i))+ggtitle(i)+theme(text = element_text(size = 24))
}
for (i in c("black","blue","yellow")) {
  ggsave(filename = paste0("A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/box/cor_",i,".pdf"),p[[i]],height = 8,width = 8)
}

### 基因与模块、表型的相关性散点图

#所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因算出相关系数， 
#如果跟性状显著相关的基因也跟某个模块显著相关，那么这些基因可能就非常重要。

#选择离散表型
nSamples<-nrow(datExpr)
levels(traitData$catgory)
choose_group <- "pro_meta" 
modNames <- substring(names(MEs), 3)
### 计算模块与基因的相关性矩阵 
## Module Membership: 模块内基因表达与模块特征值的相关性
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)
###  计算性状与基因的相关性矩阵 

## Gene significance，GS：比较样本某个基因与对应表型的相关性
## 连续型性状 直接计算
# trait <- datTraits$groupNo  
## 非连续型性状，需转为0-1矩阵, 已存于design中

trait <- as.data.frame(design[,choose_group])
geneTraitSignificance <- as.data.frame(cor(datExpr,trait,use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
names(geneTraitSignificance) <- paste0("GS")
names(GSPvalue) <- paste0("GS")
### 可视化基因与模块、表型的相关性.
selectModule <- c("black","blue","yellow")##选择自己想要的模块
pdf("A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/VSplot/step4_gene-Module-trait-significance.pdf",width=14, height=1.5*ncol(MEs))
par(mfrow=c(ceiling(length(selectModule)/2),2))
for(module in selectModule){
  column <- match(module,selectModule)
  print(module)
  moduleGenes<-moduleColors==module
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for trait",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
dev.off()

###第四部分结束
###################第五部分，WGCNA可视化：TOMplot 、Eigengene-adjacency-heatmap#####################

#[这是所谓的重新计算一遍TOM矩阵]
#dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = power)
#这是直接导入,一样的
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM

###全基因全大图
geneTree<-net$dendrograms[[1]]
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
#这一部分特别耗时，行列同时做层级聚类
TOMplot(plotTOM, geneTree, moduleColors,
        main = "Network heatmap plot, all genes")

###限制基因的数量来加快绘图
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
nSelect = 0.1*nGenes
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method = "average")
moduleColors = labels2colors(net$colors)
selectColors = moduleColors[select];
sizeGrWindow(9,9)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
#改变热图的深色背景为白色背景：
library(gplots)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
pdf(file = "A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/tomplot_heat.pdf",width = 5,height = 5)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes", col=myheatcol)
dev.off()

###绘制模块之间相关性图,Eigengene-adjacency-heatmap
MEs = net$MEs# module eigengene 可以绘制线图，作为每个模块的基因表达趋势的展示
library("stringr")
colnames(MEs) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MET = orderMEs(MEs)
pdf("A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/module_cor.pdf",width = 10,height = 10)
plotEigengeneNetworks(MET, setLabels="", 
                      marDendro = c(0,4,1,4),  # 留白：下右上左
                      marHeatmap = c(5,5,1,2), # 留白：下右上左
                      cex.lab = 0.8,
                      xLabelsAngle = 90)
dev.off()
#添加表型数据
##连续型性状，MET = orderMEs(cbind(MEs,traitData$age))
#离散型性状
design
pro_meta<-as.data.frame(design[,"pro_meta"])
names(pro_meta)<-"pro_meta"
MET = orderMEs(cbind(MEs, pro_meta))
pdf(file = "A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/prometa_cor_gene_heat.pdf",width = 8,height = 8)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", 
                      marDendro = c(1,3,2,4),
                      marHeatmap = c(5,4,1,2), 
                      plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
#要想拆分聚类图和热图，可以用以下代码实现。
sizeGrWindow(6,6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
##########第五部分结束############


########################第六部分，对感兴趣模块的基因进行GO分析#################
library("GO.db")
library("org.Hs.eg.db")
library("anRichment")
library(clusterProfiler)
#
OrgDb<-"org.Hs.eg.db"
genetype = "SYMBOL"
intModules = unique(moduleColors)
choose_module <-c("green","black","red","brown","yellow","turquoise","blue","grey")
gene_module <- data.frame(gene=colnames(datExpr),module=moduleColors)
write.csv(gene_module,file = "./gene_moduleColors.csv",row.names = F, quote = F)
tmp <- bitr(gene_module$gene,fromType = genetype,  # "SYMBOL"   "ENSEMBL"
            toType = "ENTREZID",
            OrgDb = OrgDb )
gene_module_entrz <- merge(tmp,gene_module, by.x=genetype, by.y="gene")
choose_gene_module_entrz <- gene_module_entrz[gene_module_entrz$module %in% choose_module,]
###run GO
formula_res <- compareCluster(
  ENTREZID~module,
  data = choose_gene_module_entrz,
  fun = "enrichGO",
  OrgDb = OrgDb,
  ont = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
###精简GO富集的结果,去冗余
lineage1_ego <- simplify( 
  formula_res,
  cutoff=0.5,
  by="p.adjust",
  select_fun=min
)
getwd()
save(gene_module, formula_res, lineage1_ego, file="module_GO_term.Rdata")
write.csv(lineage1_ego@compareClusterResult,file="module_GO_term.csv")
### 绘制dotplot图
load("A://HCC/Bulkseq/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/module_GO_term.Rdata")
selectModule <- c("black","blue")
dotp <- dotplot(lineage1_ego[lineage1_ego@compareClusterResult$module%in%selectModule, asis=T],
                showCategory=25,
                includeAll = TRUE, #将有overlap的结果也展示出来
                label_format=90)+theme_classic()
ggsave(filename = "A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/GO_KEGG/selectmodules_black_blue.pdf",plot = dotp,width = 8,height = 12)
#自带的GO包
library(ggplot2)
library("dplyr")
#
#source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R")
#installAnRichment()
#
GOcollection = buildGOcollection(organism = "human")
entrez = convert2entrez(organism = "human", symbol = colnames(datExpr))
table(!is.na(entrez))
GOenrichment = enrichmentAnalysis(
  classLabels = moduleColors, identifiers = entrez,
  refCollection = GOcollection,
  useBackground = "given",
  threshold = 0.05,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey")

goinput<-subset(GOenrichment$enrichmentTable,class%in%c("red","brown"))
goinput$inGroups<-factor(goinput$inGroups,levels = c("GO|GO.BP|GO","GO|GO.CC|GO","GO|GO.MF|GO"))
goinput<-arrange(goinput,inGroups,Bonferroni)
goinput<-subset(goinput,Bonferroni<0.05&FDR<0.05&inGroups=="GO|GO.BP|GO")[1:50,]
goinput$dataSetName<-factor(goinput$dataSetName,levels = rev(unique(goinput$dataSetName)))
x=goinput$class
y=goinput$dataSetName
p = ggplot(goinput,aes(x,y))
p1 = p + geom_point(aes(size=enrichmentRatio,color=-log(pValue,base = 10),shape=inGroups))+
  scale_color_gradient(low = "blue", high = "red")+scale_size(range = c(3, 10))
p2 = p1 + labs(color=expression(-log[10](pvalue)),
               size="enrichmentRatio",
               x="Module",
               y="Go_term",
               title="Go enrichment of Module Genes")+theme_bw()+theme(text = element_text(size = 20))
#######第六部分结束#######



##################第七部分，感兴趣基因模块绘制热图###############
moduleColors <- labels2colors(net$colors)
selectModule <- c("black","blue","yellow")#
for (module in selectModule) {
  color_select<-function(class_){
  color<-as.data.frame(readxl::read_xlsx("A://HCC/配色.xlsx",sheet = "Sheet1"))[,1:2]
  colnames(color)<-c("clusters","color")  
  loc<-match(class_,color$clusters)
  target<-color[loc,"color"]
  names(target)<-color[loc,"clusters"]
  return(target)
}
dat=datExpr[,moduleColors==module]
library(pheatmap)
n=t(scale(dat))
group_list=traitData$catgory
ac=data.frame(g=group_list)
rownames(ac)=colnames(n)
p<-pheatmap::pheatmap(n,
                   fontsize = 8,
                   show_colnames =T,
                   show_rownames = F,
                   cluster_cols = T,
                   annotation_col =ac,
                   width = 10,
                   height = 6,
                   angle_col=45,
                   main = paste0("module_",module,"-gene heatmap"),scale = "none")
test.clust <- cbind(cluster = cutree(p$tree_row, k = 2),n)
ann_row<-as.data.frame(test.clust[,1])
colnames(ann_row)<-"status"
ann_row[ann_row$status==1,"status"]<-"GeneSet1"
ann_row[ann_row$status==2,"status"]<-"GeneSet2"
#
labs.col <- colnames(n)
labs.col[grep("TCGA",labs.col)] <- ""
#
p<-pheatmap::pheatmap(n,
                      fontsize = 8,
                      show_colnames =T,
                      show_rownames = F,
                      cluster_cols = T,
                      annotation_col =ac,
                      width = 10,
                      height = 6,
                      angle_col=45,
                      main = paste0("module_",module,"-gene heatmap"),
                      annotation_row = ann_row,
                      annotation_colors = list(g=color_select(class_ = levels(ac$g))),
                      labels_col = labs.col,scale = "none")
ggsave(filename = paste0("A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/heatmaps/hubgenes_",module,".pdf"),plot = p,height = 9,width = 15)


GS1<-rownames(ann_row)[ann_row$status=="GeneSet1"]
GS2<-rownames(ann_row)[ann_row$status=="GeneSet2"]
write.table(x = ann_row,file = paste0("A://HCC/Bulkseq/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/DEGs/",module,".txt"),sep = "\t",row.names = T,quote = F,col.names = F)

# genelist_multi<-list(GS1,GS2)
# names(genelist_multi)<-c("GeneSet1","GeneSet2")
# source("A://HCC/Bulkseq/code/dnstram_R/GO_KEGG_multi.R")
# #
# save(list = c("result_GO_multi","result_GO_sim_multi","result_kegg_multi","p3"),file = paste0("A://HCC/Bulkseq/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/GO_KEGG/",module,".RData"))
# #
# ggsave(filename = paste0("A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/GO_KEGG/compares/hubgenes_GO_",module,".pdf"),plot = p3,height = 10,width = 15)
# #
# ggsave(filename = paste0("A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/GO_KEGG/compares/hubgenes_KEGG_",module,".pdf"),plot = p4,height = 10,width = 15)
}
load("A://HCC/Bulkseq/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/GO_KEGG/blue.RData")
library("ggplot2")
ggsave(filename = paste0("A://HCC/figs/bulk/WGCNA/HCC_mix/MADTOP5000_onlyT_h2/GO_KEGG/compares/hubgenes_GO_","blue",".pdf"),plot = p3,height = 15,width = 15)


###########第七部分，导出网络用于Cytoscape################
module = "black"
getfile<-function(module){
### 提取感兴趣模块基因名
gene <- colnames(datExpr) 
inModule <- moduleColors==module
modgene <- gene[inModule]
### 模块对应的基因关系矩阵
modTOM <- TOM[inModule,inModule]
dimnames(modTOM) <- list(modgene,modgene)
### 筛选连接度最大的top100基因
nTop = 100
IMConn = softConnectivity(datExpr[, modgene]) #计算连接度
top = (rank(-IMConn) <= nTop) #选取连接度最大的top100
filter_modTOM <- modTOM[top, top]
# for cytoscape
cyt <- exportNetworkToCytoscape(filter_modTOM,
                                edgeFile = paste("edges-", paste(module, collapse="-"), ".txt", sep=""),
                                nodeFile = paste("nodes-", paste(module, collapse="-"), ".txt", sep=""),
                                weighted = TRUE,
                                threshold = 0.1,  #weighted权重筛选阈值，可调整
                                nodeNames = modgene[top], 
                                nodeAttr = moduleColors[inModule][top])
write.table(cyt$edgeData,paste("./forcyto/edge_",module,"_100.txt"),sep = "\t",col.names = T,row.names = F,quote = F)
write.table(cyt$nodeData,paste("./forcyto/node_",module,"_100.txt"),sep = "\t",col.names = T,row.names = F,quote = F)
}
getfile(module = "green")

########################分步法，展示每一步干了啥,仅供参考#########################
### 计算邻接矩阵
adjacency = adjacency(dataExpr, power = power)

### 把邻接矩阵转换为拓扑重叠矩阵，以降低噪音和假相关，获得距离矩阵。
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

### 层级聚类计算基因之间的距离树 
geneTree = hclust(as.dist(dissTOM), method = "average")

### 模块合并
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)

### 通过计算模块的代表性模式和模块之间的定量相似性评估，合并表达图谱相似
#的模块
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
MEDissThres = 0.25

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged
## 分步法完结