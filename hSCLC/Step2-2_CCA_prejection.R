library(ggplot2)
library(Seurat)
library(dplyr)
library(stringr)
library(sctransform)
library(RColorBrewer)
library(homologene)
library(harmony)
options(future.globals.maxSize = 5 * 1024^3)

# Prepare mouse reference data
load("NOD.RData")
load("NOD_metadata.RData")

set.seed(123)
ref_anno <- metadata %>% group_by(Cluster) %>% sample_n(1000) %>% as.data.frame()
ref <- subset(Origin,cells = ref_anno$Cell)
rm(Origin,metadata);gc()

H2M <- mouse2human(genes = rownames(ref))
H2M <- H2M[which(!duplicated(H2M$mouseGene)),]
H2M <- H2M[which(!duplicated(H2M$humanGene)),]

expr <- ref@assays$SCT@counts[H2M$mouseGene,]
rownames(expr) <- H2M$humanGene
new_assay <- CreateAssayObject(counts = expr)
ref[["Anno"]] <- new_assay

ref <- SCTransform(ref,assay = "Anno",vars.to.regress = "percent.mt")
ref <- RunPCA(ref,assay = "SCT")
ref <- RunHarmony(ref,group.by.vars = "mouse",ncores = 10)
ref <- RunUMAP(ref,dims = 1:30,reduction = "harmony",reduction.name = "umap.harmony")
ref$Cluster <- Idents(ref)
ref$Cluster_merge <- str_extract(ref$Cluster,"(C\\d:\\s)(.*)(\\s\\d$)",group = 2)
save(ref,ref_anno,file = "NOD_EN_ref.RData")

load("SCLC_A.RData")
anchors <- FindTransferAnchors(reference = ref, query = SCLCa, dims = 1:30,
                               normalization.method = "SCT",reference.assay = "SCT",
                               reduction = "cca",verbose = T)
predictions <- TransferData(anchorset = anchors, refdata = ref$Cluster_merge, dims = 1:30,weight.reduction = "cca")
saveRDS(predictions,file = "Transfer_prediction_cluster_merge.rds")
SCLCa$Cluster_merge <- prediction$predicted.id
SCLCa <- SetIdent(SCLCa,value = "Cluster_merge")
SCLCa <- RunUMAP(SCLCa,dims = 1:20,reduction = "harmony",min.dist = 0.5,spread = 0.3,
                 reduction.name = "umap.harmony.new",assay = "SCT")

mycol4 <- c("Epithelial-like"="#FF6A6A","Neuronal-like"="#4DBBD5FF","Hybrid state"="#FF8C00")
pdf("TransferData_prediction_cluster_merge_score.pdf",height = 4,width = 4)
ggplot(predictions,aes(`predicted.id`,`prediction.score.max`))+
  geom_violin(aes(fill=`predicted.id`),scale = "width",show.legend = F)+
  geom_boxplot(width=0.1,outliers = F)+scale_fill_manual(values = mycol4)+theme_classic()+
  scale_y_continuous(limits = c(0,1),expand = c(0,0))
dev.off()

pdf("SCLC_A_UMAP_Cluster_merge_new_min.dist0.5_spread0.3_new.pdf",height = 5,width = 6)
DimPlot(SCLCa,reduction = "umap.harmony.new",raster = F,cols = mycol4)
dev.off()

EN_markers <- c("CDH1","CLDN7","CLDN3","KRT8","KRT18","EFNA5","NCAM1","CDH2","DPP10","ZEB1")
pdf("EN_marker_expression.pdf",height = 4,width = 6)
DotPlot(SCLCa,features = EN_markers,group.by = "Cluster_abbr")+
  scale_y_discrete(limits=rev)+scale_color_gradientn(colours = brewer.pal(7,"OrRd"))+
  scale_size_continuous(range = c(0,5),limits = c(0,100))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

save(SCLCa,mycol4,file = "SCLCa_cca.RData")