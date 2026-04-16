library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(RColorBrewer)

Origin <- readRDS("step5.integrate_splitByMouse_nod_rmImmuneEndothelial_splitByMouse_harmony_1_20_0.2.rds")
Origin$mouse <- str_remove(Origin$mouse,"C")
Origin$sampleID <- str_replace_all(Origin$sampleID,c("C4007_LLT-\\d"="4007_PT","C4007_RLT"="4007_LuM","C4007_ThyT"="4007_ThyM","C4007_LvT"="4007_LvM","C4007_OvaT"="4007_OvaM","C4007_IntT"="4007_IntM",
                                                     "C4011_LLT-\\d"="4011_PT","C4011_RLT"="4011_LuM","C4011_ThyT"="4011_ThyM","C4011_LvT"="4011_LvM","C4011_LOvaT"="4011_OvaM-1","C4011_ROvaT"="4011_OvaM-2","C4011_IntT"="4011_IntM","C4011_SpT"="4011_SpnM",
                                                     "C4033_LLT-\\d"="4033_PT","C4033_RLT"="4033_LuM","C4033_ThyT"="4033_ThyM","C4033_LvT"="4033_LvM","C4033_LVs"="4033_LvmM","C4033_Spt"="4033_SpnM"))
Origin$Organ <- str_extract(Origin$sampleID,"(?<=_).*$")
Origin$organ <- str_extract(Origin$sampleID,"(?<=_)[^-]*")
Origin$organ <- str_replace(Origin$organ,"LvmM","LvM")

mycol1 <- c('PT'="#D32F2F", 'LuM'="#1976D2", 'ThyM'="#FBC02D", 'LvM-1'="#0288D1", 'LvM-2'="#388E3C", 
            'LvM-3'="#8C8C8C", 'LvM-4'="#8E24AA",'LvM' ="#FF5722", 'LvmM'="#FB9A99", 
            'OvaM'="#8BC34A",'OvaM-1'="#A65628", 'OvaM-2'="#F44336", 'IntM'="#00BCD4", 'SpnM'="#66C2A5")

mycol2 <- c('0'="#FF6A6A",'1'="#00A087FF",'2'="#4DBBD5FF",'3'="#F39B7FFF",
            '4'="#FF8C00",'5'="#00CD00",'6'="#AB82FF",'7'="#DC0000FF",'8'="#7E6148FF")

new.cluster.ids <- c("C0: Epithelial-like 1","C1: Neuronal-like 1","C2: Neuronal-like 2","C3: Epithelial-like 2","C4: Hybrid state 1",
                     "C5: Neuronal-like 3","C6: Hybrid state 2","C7: Epithelial-like 3","C8: Epithelial-like 4")
names(new.cluster.ids) <- 0:8
names(mycol2) <- new.cluster.ids
Origin <- RenameIdents(Origin, new.cluster.ids)

mutdata <- read.table("cell_numMut_raw.txt",sep = "\t",header = T)
mutdata <- aggregate(mutdata,numMut~cellBC,mean)
metadata <- read.table("step7.metadata_nod.txt",sep = "\t",header = T)
metadata$Cluster <- Idents(Origin)
metadata$Organ <- Origin$Organ
metadata$mouse <- Origin$mouse
metadata$MetStatus <- str_replace_all(metadata$MetStatus,c("Multi-metastatic"="Multi-organ-metastatic","Single-metastatic"="Single-organ-metastatic"))
metadata$MetStatus <- factor(metadata$MetStatus,levels = c("Non-metastatic","Single-organ-metastatic","Multi-organ-metastatic","Others"))
metadata$MetGroup <- str_remove_all(metadata$MetGroup,"MetGroup\\d:")
metadata$MetGroup <- str_replace_all(metadata$MetGroup,c("Primary"="PrimarySite"))
metadata$MetGroup <- factor(metadata$MetGroup,levels = c("PrimarySite-tropic","Thymus-tropic","Liver-tropic","LungMet-tropic",
                                                         "Ovary-tropic","Spine-tropic","Others"))
colnames(metadata) <- c("Cell","harmony_cluster","sampleID","organ","clone","Clone","MetStatus","MetGroup","scFitness","scPlasticity","scMetRate",
                        "CytoTRACE2_Score","CytoTRACE2_Potency","CytoTRACE2_Relative","Cluster","Organ","mouse","meanMut")
metadata <- merge(metadata,mutdata,by.x="Cell",by.y="cellBC",all.x=T)
metadata %>% mutate(Cell=factor(`Cell`,levels = Cells(Origin))) %>% arrange(`Cell`) -> metadata

MetGroup <- read.table("step8.nod_row_annotation_data0.98.txt",sep = "\t",header = T)
MetGroup$Clone <- rownames(MetGroup)
MetGroup$MetGroup <- str_remove_all(MetGroup$MetGroup,"MetGroup\\d*:")
MetGroup$MetGroup <- str_replace(MetGroup$MetGroup,"Primary","PrimarySite")
MetGroup$MetGroup <- factor(MetGroup$MetGroup,levels = c("PrimarySite-tropic","Thymus-tropic","Liver-tropic","LungMet-tropic",
                                                         "Ovary-tropic","Spine-tropic"))
MetGroup$MetStatus <- str_replace_all(MetGroup$MetStatus,c("Multi-metastatic"="Multi-organ-metastatic","Single-metastatic"="Single-organ-metastatic"))
MetGroup$MetStatus <- factor(MetGroup$MetStatus,levels = c("Non-metastatic","Single-organ-metastatic","Multi-organ-metastatic"))

metadata <- merge(metadata,MetGroup[,c("Clone","MetGroup","MetStatus")],by=c("Clone"),all.x=T)
rownames(metadata) <- metadata$Cell
metadata <- metadata[Cells(Origin),]

group <- c("#1A9993","#709AE1","#8A9197","#D2AF81","#FD7446","#197EC0")
names(group) <- names(table(MetGroup$MetGroup))

save(Origin,file = "NOD.RData")
save(metadata,mycol1,mycol2,new.cluster.ids,group,file = "NOD_metadata.RData")
save(MetGroup,group,file = "NOD_MetGroup.RData")

