library(ggplot2)
library(Seurat)
library(dplyr)
library(stringr)
library(sctransform)
library(RColorBrewer)
library(harmony)
library(clustree)
options(future.globals.maxSize = 5 * 1024^3)

# Load the SeuratObject of 39 SCLC patients scRNA-seq in 
load("F:/Lineage tracing project/13 scRNA-seq/SCLC/data/SCLC-39.RData")
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")

# QC
sce <- subset(sce, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 20)
sce$Patient <- str_extract(sce$orig.ident,"P\\d+")
sce$Patient <- factor(sce$Patient,levels = str_sort(names(table(sce$Patient)),numeric = T))
sce$Site <- str_remove(sce$orig.ident,"^P\\d+_") %>% 
  str_replace_all(c("T$"="Tumor","T1"="Tumor","B"="Blood","N"="Normal","L"="Lymph_node")) %>% 
  factor(levels = c("Normal","Tumor","Blood","Lymph_node"))

sce <- SCTransform(sce,vars.to.regress = "percent.mt",return.only.var.genes = F)
sce <- RunPCA(sce)
sce <- FindNeighbors(sce, reduction = "pca",dims = 1:16)
sce <- FindClusters(sce,resolution = c(0.1,0.2,0.3,0.5,0.8))
sce <- RunUMAP(sce,reduction = "pca",dims = 1:16)

# Integration by Harmony
sce <- RunHarmony(sce,group.by.vars = "orig.ident",ncores = 10)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:16)
sce <- FindClusters(sce,,resolution = c(0.1,0.2,0.3,0.5,0.8))
clustree(sce)
sce <- SetIdent(sce,value = "SCT_snn_res.0.2")
Idents(sce) <- factor(Idents(sce),levels = 0:17)
sce <- RunUMAP(sce, reduction = "harmony", dims = 1:16,assay = "SCT",reduction.name = "umap.harmony")

save(sce,file = "Origin.RData")

# Find markers of each cluster
sce <- PrepSCTFindMarkers(sce)
All.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(All.markers,"Result/All.markers.xls",sep = "\t",row.names = F)
top2 <- All.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

DotPlot(sce,features = unique(top2$gene),assay = "SCT")+
  scale_color_gradientn(colors = brewer.pal(7,"OrRd"))+scale_y_discrete(limits=rev)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

# Define cell identity
new.ids <- c(`0`="SCLC",`1`="SCLC",`2`="Macrophage",`3`="SCLC",`4`="Macrophage",`5`="Alveolar cell",
             `6`="T cell",`7`="Plasma cell",`8`="Basal cell",`9`="Endothelial cell",`10`="SCLC",`11`="Fibroblast",
             `12`="SCLC",`13`="Ciliated cell",`14`="SCLC",`15`="SCLC",`16`="SCLC",`17`="SCLC")
sce <- RenameIdents(sce,new.ids)
sce$Cluster <- Idents(sce)
save(sce,file = "Origin.RData")

# Extract SCLC cells
SCLC <- subset(sce,subset = Cluster=="SCLC")
rm(sce);gc()

# QC and integration
SCLC <- SCTransform(SCLC,vars.to.regress = "percent.mt")
SCLC <- RunPCA(SCLC,assay = "SCT")
SCLC <- RunHarmony(SCLC,group.by.vars = "orig.ident",ncores = 10)
SCLC <- FindNeighbors(SCLC, reduction = "harmony",dims = 1:30)
SCLC <- FindClusters(SCLC,resolution = 0.1)
SCLC <- RunUMAP(SCLC,reduction = "harmony",dims = 1:30,reduction.name = "umap.harmony")
VlnPlot(SCLC,features = c("ASCL1","NEUROD1","POU2F3","YAP1"),pt.size = 0,ncol = 2)
save(SCLC,file = "SCLC.RData")

# Extract ASCL1+ SCLC cells
SCLCa <- subset(SCLC,idents = c("0","1","2","3","4","5","6","7","9","10","12","13"))

# QC and integration
SCLCa <- RunPCA(SCLCa,assay = "SCT")
SCLCa <- RunHarmony(SCLCa,group.by.vars = "orig.ident",ncores = 10)
SCLCa <- FindNeighbors(SCLCa, reduction = "harmony",dims = 1:30)
SCLCa <- FindClusters(SCLCa,resolution = c(0.1,0.2,0.3,0.5,0.8))
SCLCa <- SetIdent(SCLCa,value = "SCT_snn_res.0.1")
SCLCa <- RunUMAP(SCLCa,reduction = "harmony",dims = 1:30,reduction.name = "umap.harmony")
DimPlot(SCLCa,reduction="umap.harmony",label = T)

save(SCLCa,file = "SCLC_A.RData")

