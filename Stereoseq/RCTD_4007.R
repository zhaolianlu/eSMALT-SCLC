library(Seurat)
library(spacexr)
library(ggplot2)

#########################################################################
## 1) Load scRNAseq data, and make reference
#########################################################################
scdata<- readRDS("/syn2/zhaolian/3.JiLab/results/2.scRNAseq/4.integrated_byTypeReClustering/step7.scdata_node_w_scores.rds")
scdata@meta.data$seurat_clusters <- as.character(scdata@meta.data$harmony_cluster)
scdata <- DietSeurat(scdata)
# scdata <- JoinLayers(scdata)
scdata

#######################
immune <- readRDS('/syn1/liangzhen/jinhua_jilab_project/result/scRNA/cellranger/immune.rds')
immune <- subset(immune,cells = rownames(immune@meta.data)[!immune@meta.data$celltype %in% c('B_cell','T_cell')])
immune@meta.data$seurat_clusters <- immune@meta.data$celltype


fibroblast <- readRDS('/syn1/liangzhen/jinhua_jilab_project/result/scRNA/cellranger/Fibroblast_BL.rds')
fibroblast@meta.data$seurat_clusters <- 'Fibroblast'

endothelial <- readRDS('/syn1/liangzhen/jinhua_jilab_project/result/scRNA/cellranger/endothelial.rds')
endothelial@meta.data$seurat_clusters <- 'Endothelial'

##-------------------------------------------------------------
# ovary <- readRDS('/syn2/zhaolian/3.JiLab/results/5.Steroseq/reference/Ovary/raw_with_age.combined.RDS')
# granulosa <- subset(ovary,cluster.names %in% c("Granulosa A","Granulosa B"))
# granulosa@meta.data$seurat_clusters <- 'Granulosa'

# theca <- subset(ovary,cluster.names %in% c("Theca"))
# theca@meta.data$seurat_clusters <- 'Theca'

# oocytes <- subset(ovary,cluster.names %in% c("Oocytes"))
# oocytes@meta.data$seurat_clusters <- 'Oocytes'
##-------------------------------------------------------------

# scdata.merged <- merge(x=scdata,y=list(immune,fibroblast,endothelial,granulosa,theca,oocytes),
#                        add.cell.ids = c('tumor','immune','fibroblast','endothelial','granulosa','theca','oocytes'))
scdata.merged <- merge(x=scdata,y=list(immune,fibroblast,endothelial),
                       add.cell.ids = c('tumor','immune','fibroblast','endothelial'))
# scdata.merged <- JoinLayers(scdata.merged)
scdata.merged[["RNA"]] <- JoinLayers(scdata.merged[["RNA"]])

#######################
Idents(scdata.merged) <- "seurat_clusters"

set.seed(123)                      # 为了结果可复现，先固定随机种子
scdata.down <- subset(
  scdata.merged,
  downsample = 2000               # 每个 seurat_clusters 至多保留 2000 个细胞
)
scdata.merged <- scdata.down
#######################
# Idents(scdata.merged) <- scdata.merged@meta.data$seurat_clusters
# extract information to pass to the RCTD Reference function
counts <- scdata.merged[["RNA"]]$counts
cluster <- as.factor(scdata.merged@meta.data$seurat_clusters)
names(cluster) <- colnames(scdata.merged)
nUMI <- scdata.merged$nCount_RNA
names(nUMI) <- colnames(scdata.merged)
reference <- Reference(counts, cluster,nUMI)
#########################################################################
## 2) Load spatial data
#########################################################################
spatial <- readRDS("/syn2/zhaolian/3.JiLab/results/5.Steroseq/02.gefToGem/step1.st_cellbin_seurat_4007.rds")
# set up query with the RCTD function SpatialRNA
counts <- spatial[["RNA"]]$counts
coords <- spatial@meta.data[,c("x", "y")]
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

#########################################################################
## 3) run RCTD
#########################################################################
RCTD <- create.RCTD(query, reference, max_cores = 20)
RCTD <- run.RCTD(RCTD, doublet_mode = 'doublet')
saveRDS(RCTD,'/syn2/zhaolian/3.JiLab/results/5.Steroseq/03.cellbin/step2.RCTD_results_doublet_4007.rds')

RCTD <- create.RCTD(query, reference, max_cores = 20)
RCTD <- run.RCTD(RCTD, doublet_mode = 'full')
saveRDS(RCTD,'/syn2/zhaolian/3.JiLab/results/5.Steroseq/03.cellbin/step2.RCTD_results_full_4007.rds')
