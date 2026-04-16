library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(ggpubr)

load("SCLCa_cca.RData")

## Calculating CytoTRACE score
library(CytoTRACE)
expr <- as.matrix(SCLCa@assays$SCT@counts)
cells <- ncol(expr)
results <- CytoTRACE(mat = expr, subsamplesize = 1000,enableFast = T, ncores = 20)

SCLCa$CytoTRACE_score <- results$CytoTRACE
SCLCa$CytoTRACE_rank <- results$CytoTRACErank

pdf("SCLC_A_UMAP_CytoTRACE_score.pdf",height = 5,width = 6)
FeaturePlot(SCLCa,reduction = "umap.harmony.new",features = "CytoTRACE_score",raster = F)+
  scale_color_gradientn(colors = brewer.pal(9,"OrRd"))+coord_fixed()
dev.off()

pdf("SCLCa_Violinplot_CytoTRACE_score.pdf",height = ,width = 6)
ggplot(SCLCa@meta.data,aes(`Cluster_merge`,`CytoTRACE_score`))+
  geom_violin(scale="width",aes(fill=`Cluster_merge`),show.legend = F)+
  geom_boxplot(width=0.2,outliers = F)+scale_fill_manual(values = mycol4)+theme_classic()
dev.off()

## Calculating N-E score
# read the mouse reference data of N-E signature genes
NE_markers50 <- readxl::read_xlsx("NE_markers_ref_top50_top25.xlsx",sheet = 1,col_names = T) %>% as.data.frame()
expr <- SCLCa@assays$SCT@data %>% as.matrix()
NE_score <- c()
pb <- txtProgressBar(min = 1,max = ncol(expr),style = 3)
for (i in 1:ncol(expr)) {
  tmp <- expr[,i]
  n_cor <- cor.test(tmp[NE_markers50$humanGene],NE_markers50$N_expression,method = "spearman")$estimate
  e_cor <- cor.test(tmp[NE_markers50$humanGene],NE_markers50$E_expression,method = "spearman")$estimate
  NE_score[i] <- (n_cor-e_cor)/2
  setTxtProgressBar(pb,i)
  rm(tmp,n_cor,e_cor)
}
SCLCa$NE_spearman_markers50 <- NE_score

pdf("NE_spearman_markers50_by_cluster_merge.pdf",height = 4,width = 4)
ggplot(SCLCa@meta.data,aes(`Cluster_merge`,`NE_spearman_markers50`))+
  geom_violin(aes(fill=`Cluster_merge`),scale="width",show.legend = F)+
  geom_boxplot(width=0.2,outliers = F)+theme_classic()+
  scale_fill_manual(values = mycol4)
dev.off()

pdf("NE_spearman_markers50_by_Site.pdf",height = 4,width = 4)
ggplot(SCLCa@meta.data,aes(`Site`,`NE_spearman_markers50`))+
  geom_violin(aes(fill=`Site`),scale="width",show.legend = F)+
  geom_boxplot(width=0.2,outliers = F)+theme_classic()+
  scale_fill_manual(values = mycol4)
dev.off()

## Calculating Metastatic score
library(UCell)
library(homologene)
# read the Ascl+ SCLC metastatic geneset from the supplmentary table of Feifei Na et al. Nat cancer 2022
filepath <- "43018_2022_361_MOESM3_ESM.xlsx"
Ascl1_Metastatic_geneset <- readxl::read_xlsx(filepath,sheet = 6,col_names = T,skip = 1) %>% as.data.frame()
Ascl1_Metastatic_geneset <- Ascl1_Metastatic_geneset$gene_name[which(Ascl1_Metastatic_geneset$logFC>=0.25 & Ascl1_Metastatic_geneset$p_adj<=0.05)]
SCLCa <- AddModuleScore_UCell(SCLCa,features = list(ASCL1_SCLC_Metastasis=mouse2human(Ascl1_Metastatic_geneset$humanGene),assay = "SCT",ncores = 10,slot = "data")

pdf("Metastatic_score_by_cluster_merge.pdf",height = 4,width = 4)
ggplot(SCLCa@meta.data,aes(`Cluster_merge`,`ASCL1_SCLC_Metastasis_UCell`))+
  geom_violin(aes(fill=`Cluster_merge`),scale="width",show.legend = F)+
  geom_boxplot(width=0.2,outliers = F)+theme_classic()+
  scale_fill_manual(values = mycol4)
dev.off()

## Calculating the regulon score
library(pheatmap)
# read the mouse regulon information from the pySCENIC result
regulon <- read.csv("pySCENIC/reg.csv",sep = ",",header = T)
regulon <- regulon[3:nrow(regulon),]
colnames(regulon) <- c("TF","MotifID","AUC","NES","MotifSimilarityQvalue","OrthologousIdentity",
                       "Annotation","Context","TargetGenes","RankAtMax")
Top10_regulon <- read.table("E:/SMALT_scRNA-seq/Result_new_20251230/result3.rss_cellType_filtered.csv",sep = ",",header = T)
Top10_regulon <- reshape2::melt(Top10_regulon)
colnames(Top10_regulon) <- c("Cluster","Regulon","RSS")
Top10_regulon$Regulon <- str_remove_all(Top10_regulon$Regulon,"\\.*$")

# extract the TFs-target network
ppi <- list()
for (i in 1:nrow(regulon)) {
  Target <- str_remove_all(regulon$TargetGenes[i],"\\[|\\]")
  Target <- str_extract_all(Target,"\\(\\'\\w+\\',\\s[\\d\\.]+\\)") %>% unlist()
  Target <- str_remove_all(Target,"^\\(|\\)$|\\s")
  Target <- str_split(Target,",",simplify = T)
  tmp <- data.frame(TF=regulon$TF[i],Target=str_remove_all(Target[,1],"\\'"),weight=as.numeric(Target[,2]))
  ppi[[i]] <- tmp
  rm(tmp,Target)
}
ppi <- data.table::rbindlist(ppi) %>% as.data.frame()
ppi <- aggregate(ppi,weight~TF+Target,max)
M2H <- mouse2human(ppi$Target)
ppi <- merge(ppi,M2H[,c("mouseGene","humanGene")],by.x="Target",by.y="mouseGene")
regulon <- split(ppi$humanGene,ppi$TF)
regulon <- regulon[which(names(regulon) %in% Top10_regulon$Regulon)]
M2H <- mouse2human(names(regulon))
all(names(regulon)==M2H$mouseGene)
names(regulon) <- paste0(M2H$humanGene,"_target")

# Find the target genes of Zeb1
Zeb1_target <- ppi$Target[which(ppi$TF=="Zeb1")]
SCLCa <- AddModuleScore_UCell(SCLCa,features = regulon,assay = "SCT",ncores = 10,slot = "data",maxRank = 3000)

pdf("ZEB1_regulon_UCell_expression.pdf",height = 4,width = 4)
ggplot(SCLCa@meta.data,aes(`Cluster_merge`,`ZEB1_target_UCell`))+
  geom_violin(aes(fill=`Cluster_merge`),scale="width",show.legend = F)+
  scale_fill_manual(values = mycol4)+theme_classic()+
  stat_compare_means(comparisons = list(c("Epithelial-like","Neuronal-like")),method = "wilcox.test")
dev.off()

plotdata <- data.frame(SCLCa@meta.data[,paste0(names(regulon),"_UCell")])
plotdata$Cluster_merge <- SCLCa$Cluster_merge
plotdata <- reshape2::melt(plotdata,id.vars="Cluster_merge",variable.name="regulon",value.name="UCell")
plotdata2 <- aggregate(plotdata$UCell,list(Cluster=plotdata$Cluster_merge,regulon=plotdata$regulon),mean)
plotdata2 <- reshape2::dcast(plotdata2,Cluster~regulon,value.var = "x")
rownames(plotdata2) <- plotdata2$Cluster
plotdata2 <- as.matrix(plotdata2[,-1])
colnames(plotdata2) <- str_remove(colnames(plotdata2),"_target_UCell")

mycol <- colorRampPalette(brewer.pal(9,"OrRd"))(256)
pdf("SCLC_A_Regulon_UCell_heatmap.pdf",height = 4,width = 8)
pheatmap(plotdata2,scale = "column",color = mycol,border_color = NA,angle_col = 45,clustering_method = "ward.D")
dev.off()

pdf("SCLC_A_TF_expression_dotplot.pdf",height = 4,width = 10)
DotPlot(SCLCa,features = colnames(plotdata2),assay = "SCT",group.by = "Cluster_merge")+
  scale_size_continuous(range=c(0,5),limits = c(0,100))+scale_color_gradientn(colours = mycol)+
  scale_y_discrete(limits=rev)+theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()
