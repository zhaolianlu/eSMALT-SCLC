library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(aplot)
library(patchwork)
library(ggrepel)
library(colorspace)
library(RColorBrewer)
library(homologene)

load("NOD.RData")
load("NOD_metadata.RData")

Origin <- PrepSCTFindMarkers(Origin,assay = "SCT")
genes <- rownames(Origin)
genes <- genes[which(!str_detect(genes,"^mt-") & !str_detect(genes,"Target"))]

All.markers1 <- FindMarkers(Origin,features=genes,assay = "SCT",ident.1 = c("C1: Neuronal-like 1"),ident.2 = c("C0: Epithelial-like 1"),
                           logfc.threshold = 0,min.pct = 0.25,only.pos = F,recorrect_umi=F,test.use="MAST")
All.markers1$gene <- rownames(All.markers1)
All.markers1 <- arrange(All.markers1,desc(avg_log2FC))
All.markers1$col <- "grey50"
All.markers1$col[which(All.markers1$p_val_adj<=0.05 & All.markers1$avg_log2FC>=0.5)] <- "#D73027"
All.markers1$col[which(All.markers1$p_val_adj<=0.05 & All.markers1$avg_log2FC<= -0.5)] <- "#4575B4"
All.markers1$diff <- All.markers1$pct.1-All.markers1$pct.2

All.markers2 <- FindMarkers(Origin,features=genes,assay = "SCT",
                            ident.1 = c("C1: Neuronal-like 1","C2: Neuronal-like 2","C5: Neuronal-like 3"),
                            ident.2 = c("C0: Epithelial-like 1","C3: Epithelial-like 2","C7: Epithelial-like 3","C8: Epithelial-like 4"),
                            logfc.threshold = 0,min.pct = 0.25,only.pos = F,recorrect_umi=F,test.use="MAST")
All.markers2$gene <- rownames(All.markers2)
All.markers2 <- arrange(All.markers2,desc(avg_log2FC))
All.markers2$col <- "grey50"
All.markers2$col[which(All.markers2$p_val_adj<=0.05 & All.markers2$avg_log2FC>=0.5)] <- "#D73027"
All.markers2$col[which(All.markers2$p_val_adj<=0.05 & All.markers2$avg_log2FC<= -0.5)] <- "#4575B4"
All.markers2$diff <- All.markers2$pct.1-All.markers2$pct.2

All.markers <- merge(All.markers1[,c("gene","avg_log2FC","p_val_adj","col")],All.markers2[,c("gene","avg_log2FC","p_val_adj","col")],by="gene",all=T)
All.markers$ann <- "Nonsignif"
All.markers$ann[which(All.markers$col.x=="#D73027"&All.markers$col.y=="#D73027")] <- "upDEGs"
All.markers$ann[which(All.markers$col.x=="#4575B4"&All.markers$col.y=="#4575B4")] <- "dnDEGs"

pdf("DEGs_N_vs_E overlap C1_vs_C0.pdf",height = 4,width = 4)
ggplot(All.markers,aes(avg_log2FC.x,avg_log2FC.y))+geom_vhlines(xintercept = 0,yintercept = 0,linetype=4)+
  geom_point(aes(col=`ann`))+theme_bw()+
  scale_color_manual(values = c("Nonsignif"="grey70","upDEGs"="#D73027","dnDEGs"="#4575B4"))+
  geom_text_repel(data = plotdata[which(All.markers$ann!="Nonsignif" & All.markers$gene %in% hub_df$gene_name),],aes(label = `gene`),max.overlaps = 50)+
  xlab("Log2FC Neuronal-like vs Epithelial-like")+ylab("Log2FC C1 vs C0")
dev.off()

NE_markers <- All.markers$gene[which(All.markers$ann %in% c("upDEGs","dnDEGs"))]
E <- Origin@assays$SCT@data[NE_markers,which(Idents(Origin) %in% c("C0: Epithelial-like 1","C3: Epithelial-like 2",
                                                                   "C7: Epithelial-like 3","C8: Epithelial-like 4"))] %>% rowMeans()
N <- Origin@assays$SCT@data[NE_markers,which(Idents(Origin) %in% c("C1: Neuronal-like 1","C2: Neuronal-like 2",
                                                                   "C5: Neuronal-like 3"))] %>% rowMeans()
NE_markers <- data.frame(gene=NE_markers,E_expression = E,N_expression = N)
m2h <- mouse2human(NE_markers$gene)
NE_markers <- merge(NE_markers,m2h[,c("mouseGene","humanGene")],by.x="gene",by.y="mouseGene")
NE_markers$diff <- NE_markers$N_expression - NE_markers$E_expression
NE_markers <- arrange(NE_markers,desc(diff))

All.markers$avg_log2FC <- sqrt(1/2)*All.markers$avg_log2FC.x+sqrt(1/2)*All.markers$avg_log2FC.y
All.markers <- All.markers[which(All.markers$gene %in% NE_ref$gene),]
NE_markers50 <- All.markers[which(All.markers$ann %in% c("upDEGs","dnDEGs")),] %>% group_by(ann) %>% top_n(50,abs(`avg_log2FC`))
NE_markers50 <- NE_markers[which(NE_markers$gene %in% NE_markers50$gene),]
NE_markers25 <- All.markers[which(All.markers$ann %in% c("upDEGs","dnDEGs")),] %>% group_by(ann) %>% top_n(25,abs(`avg_log2FC`))
NE_markers25 <- NE_ref[which(NE_ref$gene %in% NE_markers25$gene),]
writexl::write_xlsx(list(NE_markers50=NE_markers50,NE_markers25=NE_markers25),"NE_markers_ref_top50_top25.xlsx")
