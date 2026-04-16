library(ggplot2)
library(stringr)
library(igraph)
library(ggraph)
library(tidygraph)
library(grid)
library(ggplotify)
library(cowplot)
library(patchwork)
library(treeio)
library(ggtree)
library(pheatmap)
library(RColorBrewer)
library(colorRamp2)

load("NOD_MetGroup.RData")
load("NOD_metadata.RData")
rownames(metadata) <- str_remove(metadata$Cell,"-\\d*$")
tree <- metadata[which(metadata$scFitness!="Others"),]
tree <- aggregate(tree,Cell~MetStatus+MetGroup+Clone+organ+Organ,length)
ann <- aggregate(tree,Cell~MetStatus+MetGroup+Clone,sum)
mycol3 <- mycol1[names(mycol3)]

filedir <- "6.fitchCount_5cells/"
files <- list.files(filedir)
cellnum <- files[which(str_detect(files,"_cellNumber\\.csv"))]
mtx <- files[which(str_detect(files,"transitionCount\\.csv"))]
phylo <- files[which(str_detect(files,"pruned\\.nwk$"))]

Origin <- list()
Origin.scaled <- list()
node_ann <- list()
Phylo <- list()
for (i in 1:length(cellnum)) {
  tmp <- read.csv(paste0(filedir,mtx[i]),
                  sep = ",",header = T,row.names = 1,check.names = F) %>% as.matrix()
  mtx0 <- matrix(0,nrow = length(mycol3),ncol = length(mycol3)) %>% `dimnames<-`(list(names(mycol3),names(mycol3)))
  mtx0[rownames(tmp),colnames(tmp)] <- tmp
  Clone <- str_extract(cellnum[i],"^\\d{4}_\\d*")
  Origin[[Clone]] <- mtx0
  
  tmp <- read.table(paste0(filedir,cellnum[i]),
                    sep = ",",header = T,check.names = F)
  colnames(tmp) <- c("Organ","Cell","Ratio")
  tmp$Organ <- factor(tmp$Organ,levels=names(mycol3))
  tmp$Clone <- Clone
  tmp$MetStatus <- ann$MetStatus[which(ann$Clone==Clone)]
  tmp$MetGroup <- ann$MetGroup[which(ann$Clone==Clone)]
  mtx1 <- mtx0
  diag(mtx0) <- 0
  tmp$WiD <- colSums(mtx0[,tmp$Organ])
  tmp$WoD <- rowSums(mtx0[tmp$Organ,])
  
  mtx1 <- mtx1/(rowSums(mtx1) %*% t(rep(1,ncol(mtx1))))
  mtx1[is.nan(mtx1)] <- 0
  Origin.scaled[[Clone]] <- mtx1
  diag(mtx1) <- 0
  tmp$WiD_rescale <- colSums(mtx1[,tmp$Organ])
  tmp$WoD_rescale <- rowSums(mtx1[tmp$Organ,])
  node_ann[[Clone]] <- tmp
  
  tmp <- read.tree(paste0(filedir,phylo[i]))
  tree_ann <- metadata[tmp$tip.label,"organ"]
  tree_ann[as_tibble(tmp)$parent] <- "undefined"
  tmp <- full_join(tmp,data.frame(node=1:(length(tmp$edge.length)+1),Organ=tree_ann),by="node")
  Phylo[[Clone]] <- tmp
  rm(tmp,Clone,mtx0,mtx1)
}

# construct a 0 matrix with nrow and ncol == clone number
mtx <- matrix(0,ncol = length(Origin.scaled),nrow = length(Origin.scaled)) %>% 
  `dimnames<-`(list(names(Origin.scaled),names(Origin.scaled)))
for (x in names(Origin.scaled)) {
  for (y in names(Origin.scaled)) {
    tmp1 <- Origin.scaled[[x]]
    diag(tmp1) <- 0
    tmp2 <- Origin.scaled[[y]] # each rescaled probability matrix has the same dimensions and row/column order
    diag(tmp2) <- 0
    mtx[x,y] <- sum(tmp1*tmp2)/(sum(tmp1^2)+sum(tmp2^2)-sum(tmp1*tmp2)) # Calculate the extended Jaccard distance
    rm(tmp1,tmp2)
  }
}
mtx <- mtx[which(!is.nan(rowSums(mtx))),which(!is.nan(colSums(mtx)))] # remove rows and columns contain NaN

col <- colorRamp2(c(0,0.3,0.6,1),c("#FCEBD5","#F99B75","#E0325E","#50193BFF"))
pattern <- c(Pattern_1="#00DCBA",Pattern_2="#F0AD00",Pattern_3="#FF80DE",Pattern_4="#00D4FF",
             Pattern_5="#7ACE00",Pattern_6="#C8A3FF",Pattern_7="#82B7FF",Pattern_8="#FF83FF",
             Pattern_9="#00C8FF",Others="#898989")
ann2 <- data.table::rbindlist(node_ann) %>% as.data.frame()
ann2 <- aggregate(ann2,Cell~Clone+MetStatus+MetGroup,sum)
rownames(ann2) <- ann2$Clone
pdf("EJaccard_distance_heatmap_ward.D_average_complete.pdf",height = 10,width = 12)
c1 <- hclust(as.dist(1-mtx),method = "ward.D")
pattern <- cutree(c1,k = 12)
ann2 <- merge(ann2,data.frame(Clone=names(pattern),Pattern=pattern),by="Clone")
rownames(ann2) <- ann2$Clone
ann2 <- ann2[c1$order,]
ann2$Pattern <- factor(paste0("pattern_",ann2$Pattern),levels = paste0("pattern_",1:12))
ann2$Pattern <- str_replace_all(as.character(ann2$Pattern),`names<-`(c(paste0("Pattern_",1:8),"Others",rep("Pattern_9",3)),paste0(unique(ann2$Pattern),"$")))
pheatmap(1-mtx,cluster_rows = c1,cluster_cols = c1,border_color = NA,annotation = ann2[,c("Pattern","MetGroup","MetStatus")],
         annotation_colors = list(MetGroup=group,MetStatus=c(`Non-metastatic`="#146C94",`Single-organ-metastatic`="#5CCCF4",
                                                             `Multi-organ-metastatic`="#F58DA9"),Pattern=pattern),
         color = colorRampPalette(rev(brewer.pal(9,"OrRd")))(256),angle_col = 45,main = paste("k=",12))
dev.off()
write.table(ann2[c1$labels[c1$order],c("Clone","MetStatus","MetGroup","Pattern")],"Cutree_group.txt",quote = F,row.names = F,sep = "\t")

# Calculate the Betweenness of each node in each clone
for (i in c1$labels[c1$order]) {
  print(i)
  tmp <- Origin[[i]]
  diag(tmp) <- 0
  tmp <- reshape2::melt(tmp)
  colnames(tmp) <- c("Source","Target","Frequency")
  tmp <- tmp[which(tmp$Frequency>0),]
  
  tmp <- merge(tmp,node_ann[[i]][,c("Organ","Cell")],by.x="Source",by.y = "Organ",all.x = T)
  tmp$group <- tmp$Source
  g <- graph_from_data_frame(tmp,vertices = node_ann[[i]])
  V(g)$label <- paste0(V(g)$name,"\n(",node_ann[[i]]$Cell,")")
  
  node_ann[[i]]$Betweenness <- betweenness(g,weights = 1/E(g)$Frequency)
  node_ann[[i]]$Normalized_Betweenness <- betweenness(g,weights = 1/E(g)$Frequency,normalized = T)
  
  tmp <- Origin.scaled[[i]]
  diag(tmp) <- 0
  tmp <- reshape2::melt(tmp)
  colnames(tmp) <- c("Source","Target","Frequency")
  tmp <- tmp[which(tmp$Frequency>0),]
  
  tmp <- merge(tmp,node_ann[[i]][,c("Organ","Cell")],by.x="Source",by.y = "Organ",all.x = T)
  tmp$group <- tmp$Source
  g <- graph_from_data_frame(tmp,vertices = node_ann[[i]])
  V(g)$label <- paste0(V(g)$name,"\n(",node_ann[[i]]$Cell,")")
  
  node_ann[[i]]$Betweenness_rescaled <- betweenness(g,weights = log10(1/E(g)$Frequency))
  node_ann[[i]]$Normalized_Betweenness_rescaled <- betweenness(g,weights = log10(1/E(g)$Frequency),normalized = T)
  
  rm(g,tmp)
}

# Integrate the metastatic probability matrix of each mouse
summ <- list()
summ.scaled <- list()
for (i in c("4007","4011","4033")) {
  clone <- names(Origin)[which(str_detect(names(Origin),i))]
  mtx0 <- matrix(0,ncol = length(mycol3),nrow = length(mycol3)) %>% `dimnames<-` (list(names(mycol3),names(mycol3)))
  for (j in clone) {mtx0 <- mtx0+Origin[[j]]}
  summ[[i]] <- mtx0
  mtx0 <- mtx0/(rowSums(mtx0) %*% t(rep(1,length(mycol3))))
  mtx0[is.nan(mtx0)] <- 0
  summ.scaled[[i]] <- mtx0
  rm(mtx0,clone)
}

plotdata <- data.table::rbindlist(node_ann) %>% as.data.frame()
plotdata$Mouse <- str_extract(plotdata$Clone,"^\\d{4}")

pdf("Weighted degree MetGroup.pdf",height = 3,width = 6)
ggplot(plotdata,aes(Organ,WiD_rescale))+geom_boxplot(aes(fill=MetGroup),show.legend = F,outliers = F)+
  geom_jitter(width = 0.1,shape=1)+facet_grid(.~`MetGroup`,scale="free_x",space="free_x")+theme_classic()+
  scale_fill_manual(values = group)+ylab("In-strenghth")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),strip.background = element_blank())
ggplot(plotdata,aes(Organ,WoD_rescale))+geom_boxplot(aes(fill=MetGroup),show.legend = F,outliers = F)+
  geom_jitter(width = 0.1,shape=1)+facet_grid(.~`MetGroup`,scale="free_x",space="free_x")+theme_classic()+
  scale_fill_manual(values = group)+ylab("Out-strenghth")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),strip.background = element_blank())
dev.off()

pdf("Weighted degree.pdf",height = 3,width = 3)
ggplot(plotdata,aes(Organ,WiD_rescale))+geom_boxplot(aes(fill=Organ),show.legend = F,outliers = F)+
  geom_jitter(width = 0.1,shape=1)+theme_classic()+
  scale_fill_manual(values = mycol3)+ylab("In-strenghth")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggplot(plotdata,aes(Organ,WoD_rescale))+geom_boxplot(aes(fill=Organ),show.legend = F,outliers = F)+
  geom_jitter(width = 0.1,shape=1)+theme_classic()+
  scale_fill_manual(values = mycol3)+ylab("Out-strenghth")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),strip.background = element_blank())
dev.off()

# Calculate the Betweenness of integrated metastatic probability matrices
node_ann_summ <- aggregate(plotdata,Cell~Mouse+Organ,sum)
node_ann_summ <- split(node_ann_summ,node_ann_summ$Mouse)
for (i in names(summ)) {
  tmp <- summ[[i]]
  diag(tmp) <- 0
  tmp <- reshape2::melt(tmp)
  colnames(tmp) <- c("Source","Target","Frequency")
  tmp <- tmp[which(tmp$Frequency!=0),]
  tmp$group <- tmp$Source
  write.table(tmp,paste0("mtx/",i,"_integrated_count.txt"),sep = "\t",quote = F,row.names = F)
  
  g <- graph_from_data_frame(tmp,vertices = node_ann_summ[[i]][,c("Organ","Cell")])
  V(g)$label <- paste0(V(g)$name,"(",V(g)$Cell,")")
  
  node_ann_summ[[i]] <- merge(node_ann_summ[[i]],aggregate(tmp,Frequency~Target,sum),by.x="Organ",by.y="Target",all.x=T)
  node_ann_summ[[i]] <- merge(node_ann_summ[[i]],aggregate(tmp,Frequency~Source,sum),by.x="Organ",by.y="Source",all.x=T)
  colnames(node_ann_summ[[i]])[c(4,5)] <- c("WiD","WoD")
  Betweenness <- data.frame(Betweenness=betweenness(g,directed = T,weights =1/tmp$Frequency),
                            Normalized_Betweenness=betweenness(g,directed = T,weights = 1/tmp$Frequency,normalized = T),
                            Organ=V(g)$name)
  node_ann_summ[[i]] <- merge(node_ann_summ[[i]],Betweenness,by="Organ",all.x = T)
  
  tmp <- summ.scaled[[i]]
  diag(tmp) <- 0
  tmp <- reshape2::melt(tmp)
  colnames(tmp) <- c("Source","Target","Frequency")
  tmp <- tmp[which(tmp$Frequency!=0),]
  tmp$group <- tmp$Source
  write.table(tmp,paste0("mtx/",i,"_integrated_Probability.txt"),sep = "\t",quote = F,row.names = F)
  
  g <- graph_from_data_frame(tmp,vertices = node_ann_summ[[i]][,c("Organ","Cell")])
  V(g)$label <- paste0(V(g)$name,"(",V(g)$Cell,")")
  
  node_ann_summ[[i]] <- merge(node_ann_summ[[i]],aggregate(tmp,Frequency~Target,sum),by.x="Organ",by.y="Target",all.x=T)
  node_ann_summ[[i]] <- merge(node_ann_summ[[i]],aggregate(tmp,Frequency~Source,sum),by.x="Organ",by.y="Source",all.x=T)
  colnames(node_ann_summ[[i]])[c(8,9)] <- c("WiD_rescaled","WoD_rescaled")
  Betweenness <- data.frame(Betweenness_rescaled=betweenness(g,directed = T,weights =log10(1/tmp$Frequency)),
                            Normalized_Betweenness_rescaled=betweenness(g,directed = T,weights = log10(1/tmp$Frequency),normalized = T),
                            Organ=V(g)$name)
  node_ann_summ[[i]] <- merge(node_ann_summ[[i]],Betweenness,by="Organ",all.x = T)
  
  tmp <- node_ann_summ[[i]][,c("Organ","Mouse","Cell")]
  tmp$label <- paste0(tmp$Organ,"(",tmp$Cell,")")
  write.table(tmp,paste0("mtx/",i,"_node_annotation.txt"),sep = "\t",quote = F,row.names = F)
  
  rm(g,tmp,Betweenness)
}

# Output the clonal metastatic network information
for (i in names(node_ann)) {
  tmp <- Origin.scaled[[i]]
  tmp <- tmp[which(rowSums(tmp)!=0),which(colSums(tmp)!=0),drop=F]
  tmp <- reshape2::melt(tmp)
  colnames(tmp) <- c("Source","Target","Probability")

  tmp <- tmp[which(tmp$Probability>0),]
  tmp$from <- tmp$Source
  tmp$label <- round(tmp$Probability,3)
  write.table(tmp,file = paste0("mtx/",i,"_Transition_Probability.txt"),sep = "\t",quote = F,row.names = F)
  
  tmp <- node_ann[[i]]
  tmp <- tmp[,c(1,2,3,4)]
  tmp$label <- paste0(tmp$Organ,"(",tmp$Cell,")")
  write.table(tmp,file = paste0("mtx/",i,"_node_annotation.txt"),sep = "\t",quote = F,row.names = F)
  rm(tmp)
}

plotdata2 <- list()
for (i in c1$labels[c1$order]) {
  tmp <- Origin.scaled[[i]]
  diag(tmp) <- 0
  tmp <- reshape2::melt(tmp)
  colnames(tmp) <- c("Source","Target","Probability")
  tmp$Clone <- i
  plotdata2[[i]] <- tmp
}
plotdata2 <- data.table::rbindlist(plotdata2) %>% as.data.frame()
plotdata2 <- plotdata2[which(plotdata2$Source!=plotdata2$Target),]
plotdata2 <- merge(plotdata2,ann2[,c("Clone","Pattern")],by="Clone")
plotdata3 <- aggregate(plotdata2,Probability~Pattern+Source+Target,mean)
plotdata3_summ <- aggregate(plotdata2,Probability~Pattern+Source+Target,function(x){length(which(x!=0))})
plotdata3 <- merge(plotdata3,plotdata3_summ,by=c("Pattern","Source","Target"))
colnames(plotdata3)[4:5] <- c("Probability","Clone")
plotdata3$Path <- paste(plotdata3$Source,plotdata3$Target,sep = " -> ")
plotdata3 <- plotdata3[which(plotdata3$Pattern!="Others"),]
Path_summ <- aggregate(plotdata3,Probability~Path,max)
plotdata3 <- plotdata3[which(plotdata3$Path %in% Path_summ$Path[which(Path_summ$Probability>=0.01)] & plotdata3$Clone!=0),]
pdf("Filter/Featured_metastasis_path_pattern.pdf",height = 4,width = 8)
ggplot(plotdata3,aes(`Path`,`Pattern`))+geom_point(aes(col=`Probability`,size=`Clone`))+
  facet_grid(.~`Source`,scale="free_x",space="free_x")+theme_bw()+
  scale_color_viridis_c()+scale_y_discrete(limits=rev)+scale_size_continuous(range = c(2,5))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),strip.background = element_blank(),strip.text = element_blank())
dev.off()

plotdata4 <- data.table::rbindlist(node_ann_summ) %>% as.data.frame()
pdf("Filter/Intgrated_degrees_mouse.pdf",height = 6,width = 4)
p1 <- ggplot(plotdata4,aes(`Organ`,`WiD_rescaled`))+geom_col(aes(fill=`Organ`),show.legend = F)+
  facet_grid(.~`Mouse`,scale="free_x",space="free_x")+scale_fill_manual(values = mycol3)+
  theme_classic()+theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))
p2 <- ggplot(plotdata4,aes(`Organ`,`WoD_rescaled`))+geom_col(aes(fill=`Organ`),show.legend = F)+
  facet_grid(.~`Mouse`,scale="free_x",space="free_x")+scale_fill_manual(values = mycol3)+
  theme_classic()+theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))
p3 <- ggplot(plotdata4,aes(`Organ`,`Betweenness_rescaled`))+geom_col(aes(fill=`Organ`),show.legend = F)+
  facet_grid(.~`Mouse`,scale="free_x",space="free_x")+scale_fill_manual(values = mycol3)+
  theme_classic()+theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,hjust = 1))
ggarrange(p1,p2,p3,ncol = 1,align = "hv")
dev.off()

save(Origin,Origin.scaled,ann,ann2,c1,node_ann,node_ann_summ,Phylo,summ,summ.scaled,plotdata,file="Route_Analysis.RData")

