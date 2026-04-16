library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(ggrepel)
library(colorspace)
library(RColorBrewer)
library(WGCNA)
library(hdWGCNA)

load("NOD.RData")
load("NOD_metadata.RData")
Origin$Cluster <- Idents(Origin)
vargene <- VariableFeatures(Origin)
vargene <- vargene[which(!str_detect(vargene,"^mt-") & !str_detect(vargene,"Target"))]

# hdWGCNA
Origin <- SetupForWGCNA(Origin,features = vargene,wgcna_name = "var3000")
Origin <- MetacellsByGroups(Origin,group.by = "Cluster",assay = "SCT",reduction = 'harmony',
                            k = 50,min_cells = 50,max_shared = 10,ident.group = "Cluster")
metacell <- GetMetacellObject(Origin)

Origin <- NormalizeMetacells(Origin)
Origin <- ScaleMetacells(Origin, features=vargene)
Origin <- RunPCAMetacells(Origin, features=vargene)
Origin <- RunUMAPMetacells(Origin,dims=1:15)
Origin <- SetDatExpr(Origin,assay = 'SCT',layer = 'data')
Origin <- TestSoftPowers(Origin,networkType = "signed")

pdf("soft threshold test for hdWGCNA.pdf",height = 6,width = 8)
p <- PlotSoftPowers(Origin)
wrap_plots(p, ncol=2)
dev.off()

Origin <- ConstructNetwork(Origin,overwrite_tom = T,soft_power = 9,deepSplit = 4,mergeCutHeight = 0.3,tom_name = "mergeCutHeight0.2")
pdf("Module_dendrogram_hdWGCNA_mergeCutHeight0.3.pdf",height = 3,width = 5)
PlotDendrogram(Origin, main='hdWGCNA Dendrogram')
dev.off()

Origin <- ModuleEigengenes(Origin,assay = "SCT")
Origin <- ModuleConnectivity(Origin,assay = "SCT",layer = "data")
Origin <- ResetModuleNames(Origin,new_name = "Module")
PlotKMEs(Origin, ncol=4)
hMEs <- GetMEs(Origin)

modules <- GetModules(Origin)
TOM <- GetTOM(Origin)
dissTOM <- 1-TOM
d <- as.dist(dissTOM)
c <- hclust(d,method = "average")
pdf("hdWGCNA_TOM_modules_heatmap_all_genes.pdf",height = 5,width = 5)
TOMplot(dissim = dissTOM^7,dendro = c,main="hdWGCNA network heatmap plot",
        Colors = modules$color,col=colorRampPalette(rev(brewer.pal(9,"OrRd")))(256))
dev.off()

Origin <- RunModuleUMAP(Origin,n_hubs = 2,n_neighbors=10,min_dist=0.1)
ModuleUMAPPlot(Origin,edge.alpha=0.25,sample_edges=T,edge_prop=0.1,label_hubs=5,keep_grey_edges=F)
hub_df <- GetHubGenes(Origin, n_hubs = 20)
write.table(modules,"hdWGCNA_modules.txt",sep = "\t",quote = F,row.names = F)
save(Origin,file = "hdWGCNA_object.RData")

# GO enrichment of hdWGCNA modules
library(clusterProfiler)
library(enrichplot)
modules <- GetModules(Origin)
modules <- subset(modules,module!="grey")
Module_genes <- split(modules$gene_name,modules$module)

GO <- list()
for(i in names(Module_genes)[-1]){
  print(i)
  input <- Module_genes[[i]]
  eg <- bitr(input,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
  go <- enrichGO(gene = unique(eg$ENTREZID),pvalueCutoff = 1,qvalueCutoff = 1,pAdjustMethod = "BH",minGSSize = 5,maxGSSize = 500,
                 OrgDb = "org.Mm.eg.db",ont = "ALL",keyType = "ENTREZID",readable = T)
  go <- as.data.frame(go)
  go$Module <- i
  GO[[i]] <- go
  rm(go)
}
GO <- data.table::rbindlist(GO) %>% as.data.frame()
GO$GeneRatio <- as.numeric(str_split(GO$GeneRatio,"\\/",simplify = T)[,1])/as.numeric(str_split(GO$GeneRatio,"\\/",simplify = T)[,2])

GO_subset <- GO[which(GO$p.adjust<=0.05 & GO$ONTOLOGY!="MF"),]

GO_top <- group_by(GO_subset,Module,ONTOLOGY) %>% top_n(-2,`p.adjust`) %>% top_n(2,`FoldEnrichment`) %>% 
  arrange() %>% as.data.frame()
ID <- GO_top$Description
GO_top <- GO_subset[which(GO_subset$Description %in% ID),]
GO_top$Module <- factor(GO_top$Module,levels = names(Module_genes)[-1])
GO_top$Description <- factor(GO_top$Description,levels = rev(unique(ID)))

pdf("Enrichment_hdWGCNA_module_genes.pdf",height = 6,width = 6)
ggplot(GO_top,aes(`Module`,`Description`))+geom_point(aes(col=`p.adjust`,size=`GeneRatio`))+theme_bw()+
  scale_color_gradientn(colors=colorRampPalette(brewer.pal(11,"RdYlBu"))(100))+
  scale_size_continuous(range = c(1,5))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

# Calculating Module UCell score
library(UCell)
Module_genes <- split(modules$gene_name,modules$module)
Origin <- AddModuleScore_UCell(Origin,features = Module_genes[-1],ncores = 10,assay = "SCT",slot = "data")

pdf("hdWGCNA_module_dotplot_cluster.pdf",height = 4,width = 8)
DotPlot(Origin,features = paste0("Module",1:10,"_UCell"),scale = T)+
  scale_color_gradientn(colors = brewer.pal(9,"OrRd"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
dev.off()

save(Origin,file = "hdWGCNA_object.RData")
write.table(Origin@meta.data[,c("Cell",paste0("Module",1:10,"_UCell"))],"hdWGCNA_module_expression.xls",
            sep = "\t",quote = F,row.names = F)

metadata <- mutate(metadata,Site=ifelse(Organ=="PT","PT","MT"))
metadata$Site <- factor(metadata$Site,levels = c("PT","MT"))
Origin@meta.data <- cbind(Origin@meta.data,metadata[,setdiff(colnames(metadata),colnames(Origin@meta.data))])
subOrigin <- subset(Origin,MetGroup!="Others")
# subOrigin <- subset(subOrigin,cells=Cells(subOrigin)[which(!(subOrigin$MetStatus=="Non-metastatic" & subOrigin$Site=="MT"))])

library(ggradar)
library(scales)
plotdata <- subOrigin@meta.data[,c("Cell","MetStatus","MetGroup","Site","Cluster",paste0("Module",1:10,"_UCell"))]
plotdata2 <- reshape2::melt(plotdata,id.vars=c("Cell","MetStatus","MetGroup","Site","Cluster"),variable.name="Module",value.name = "Module_score")
plotdata2 <- reshape2::dcast(plotdata2,MetGroup+Site~Module,value.var = "Module_score",mean)
plotdata2 <- cbind(plotdata2[,1:2],scale(plotdata2[,3:ncol(plotdata2)],center = T,scale = T))

library(fmsb)
ann <- matrix(rep(c(4,-4),10),nrow = 2,ncol = 10)
ann <- as.data.frame(ann,row.names = c("Max","Min"))
colnames(ann) <- paste0("Module",1:10)

mycol <- c("#377EB8","#E41A1C")
pdf("radarplot_hdWGCNA_only_MetGroup_fmsb.pdf",height = 5,width = 6)
for (i in names(table(plotdata2$MetGroup))) {
  tmp <- plotdata2[which(plotdata2$MetGroup==i),-1]
  rownames(tmp) <- tmp$Site
  tmp <- tmp[,-1]
  colnames(tmp) <- str_remove(colnames(tmp),"_UCell")
  tmp <- rbind(ann,tmp)
  radarchart(
    tmp, axistype = 1,
    pcol = mycol, pfcol = scales::alpha(mycol,0.3), plwd = 3, plty = 1,
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    axislabcol = "black", 
    vlcex = 1,
    caxislabels = c(-4,-2,0,2,4),title = i)
  # Add an horizontal legend
  legend(
    x = "right",legend = rownames(tmp[-c(1,2),]), horiz =F,
    bty = "n", pch = 20 , col = mycol,
    text.col = "black", cex = 1, pt.cex = 1.5
  )
}
dev.off()