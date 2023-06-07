library(fgsea)
library("ggplot2")
###GSEA
DEG_res<-read.table("A://HCC/Bulkseq/degresult/fiveclasses/nofilter_DEGs/catgory_pro_meta_vs_pro_T.txt",header = T,row.names = 1)
FC<-cbind(rownames(DEG_res),as.data.frame(DEG_res))
colnames(FC)[1]<-"SYMBOL"
FC<-dplyr::select(FC,"SYMBOL","log2FoldChange","padj")
FCthreup<- 2
FCthredn<- -2
FC$"status"<-ifelse(FC$padj>0.05,"nonsig",
                    ifelse(FC$log2FoldChange > FCthreup,"Up",ifelse(FC$log2FoldChange < FCthredn,"Dn","Stable")))
sortdf<-FC[order(FC$log2FoldChange, decreasing = T),]
gene.expr = sortdf$log2FoldChange
names(gene.expr) <- sortdf$SYMBOL

hall50<-gmtPathways(gmt.file = "A://HCC/Bulkseq/gseafile/h.all.v7.5.1.symbols.gmt")
Tumor_metastasis<-gmtPathways(gmt.file = "A://HCC/Bulkseq/gseafile/Tumor_metastasis.gmt")
allimmu<-gmtPathways(gmt.file = "A://HCC/Bulkseq/gseafile/c7.all.v2022.1.Hs.symbols.gmt")
Cell_proli<-gmtPathways(gmt.file = "A://HCC/Bulkseq/gseafile/CELL_PROLIFERATION_GO_0008283.v2022.1.Hs.gmt")
#
library("dplyr")
fgseaRes <- fgsea(pathways = hall50, 
                  stats    = gene.expr,
                  minSize  = 15,
                  maxSize  = 500,
                  eps = 0.0 )%>% as_tibble() %>% arrange(padj)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(pval) %>% 
  DT::datatable()
fgseaResTidy_sig<-subset(fgseaResTidy,padj<0.05)
#
p1<-ggplot(fgseaResTidy_sig, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES>0)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="meta_vs_T") + 
  theme_classic()+theme(text = element_text(size = 20))
#ggsave(filename = "A://HCC/figs/bulk/GSEA/GSEA_hall50.pdf",plot = p1,width = 15,height = 8)

