library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
options(future.globals.maxSize = 5 * 1024^3)
args<-commandArgs(TRUE)
my.organ=args[1]
print(my.organ)

##############################################################
## 1. Load data
##############################################################
spatial <- readRDS('/syn2/zhaolian/3.JiLab/results/5.Steroseq/03.cellbin/step5.cellbin_4007_anno.rds')
DefaultAssay(spatial) <- 'RNA'
spatial <- NormalizeData(spatial)

spatial <- subset(spatial,spot_class != 'reject' & organ== my.organ )
Idents(spatial) <- spatial@meta.data$first_typeNE
data.input = Seurat::GetAssayData(spatial, slot = "data", assay = "SCT") 
meta = data.frame(labels = Idents(spatial),row.names = names(Idents(spatial)))
spatial.locs <- as.matrix(spatial@meta.data[,c('x','y')])
message("step 1. is done...")
##############################################################
## 2. define spatial factors
##############################################################
conversion.factor = 0.5
spot.size = 16 # use the typical human cell size
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)
d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
print(min(d.spatial[d.spatial!=0]))
# spatial.factors <- list(spot.diameter = 8, spot = 10,ratio = 0.5,tol=12.5)
# spatial.factors <- list(spot.diameter = 25, spot = 10,ratio = 0.5,tol=12.5)
message("step 2. is done...")
##############################################################
## 3. create cellchat
##############################################################
cellchat <- createCellChat(object = spatial, group.by = "first_typeNE", assay = "SCT",datatype = "spatial", 
                           coordinates = spatial.locs, spatial.factors = spatial.factors)
## database
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB)
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
message("step 3. is done...")
##############################################################
## 4. run cellchat
##############################################################
future::plan("multisession", workers = 10) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = TRUE, interaction.range = 50, scale.distance = 0.5,
                              contact.dependent = TRUE, contact.range = 30)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = paste0('/syn2/zhaolian/3.JiLab/results/5.Steroseq/06.cellchat/cellchat_4007_',my.organ,'.rds'))
