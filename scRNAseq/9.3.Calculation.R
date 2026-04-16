library(jpeg)
library(EBImage)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(stringr)
library(dplyr)
library(grid)
library(ggpmisc)
library(tiff)
library(ggrepel)

# function to generate gaussian kernel
gaussian_kernel <- function(sigma,size){
  kernel <- matrix(0, nrow = size, ncol = size)
  for (x in 1:size) {
    for (y in 1:size) {
      kernel[x, y] <- exp(-((x - (size+1)/2)^2 + (y - (size+1)/2)^2) / (2 * sigma^2))
    }
  }
  kernel <- kernel/sum(kernel)
  return(kernel)
}

# Merge images and masks:
filelist <- list.files("Merge/")
filelist <- str_remove(filelist,"\\.tif")
dir.create("Segmented")
for (i in filelist) {
  print(i)
  img <- readImage(paste0("Merge/",i,".tif"))
  mask <- readTIFF(paste0("Mask/",i,"_cp_masks.tif"),as.is = T) %>% transpose()
  paintObjects(mask,img,col = "white",thick = T) %>% writeImage(paste0("Segmented/",i,"_segmented.tif"))
  rm(img,mask)
}

# Create a mask file with color labeling debris and normal cells manually, saved in Label/ directory

# Preprocess the segmented images to obtain the cell size
Result <- list()
for (i in filelist) {
  print(i)
  img <- readImage(paste0("Merge/",i,".tif"))
  mask <- readTIFF(paste0("Mask/",i,"_cp_masks.tif"),as.is = T) %>% transpose()
  tmp <- computeFeatures.shape(x = mask) %>% as.data.frame() # calculate the size information
  tmp$cell <- rownames(tmp)
  tmp$sample <- i
  
  label <- readImage(paste0("Label/",i,"_segmented.tif")) %>% EBImage::channel("grey") # read the label mask file
  label_express <- computeFeatures.basic(mask,ref = label) # calculate the label color strength within each segmented cell
  tmp$label <- ifelse(label_express[,1]<1,0,1) # add the label info to each segmented cell, with 0 as cells should be removed, 1 as cells should be kept
  tmp$annotation <- "others"
  tmp$annotation[which(tmp$s.radius.mean>=10 & tmp$s.area>=500 & tmp$label==1)] <- "tumor cells" # annotate cells with cell size and label info.
  
  Result[[i]] <- tmp
  rm(img,mask,tmp,label,label_express)
}
save(Result,file="Metadata.RData")

# Calculate the immunofluorescence intensity within each cell
for (i in filelist) {
  print(i)
  sample <- str_split(i,"_") %>% unlist()
  DAPI <- readImage(paste0("Image/",sample[1],"/",sample[2],"/",sample[2],"_ch03.tif"))
  Tubb3 <- readImage(paste0("Image/",sample[1],"/",sample[2],"/",sample[2],"_ch02.tif"))
  Cdh1 <- readImage(paste0("Image/",sample[1],"/",sample[2],"/",sample[2],"_ch01.tif"))
  Zeb1 <- readImage(paste0("Image/",sample[1],"/",sample[2],"/",sample[2],"_ch00.tif"))
  mask <- readTIFF(paste0("Mask/",i,"_cp_masks.tif"),as.is = T) %>% transpose()
  
  # calculate the coordinates and staining intensity within each segmented cell
  loc <- computeFeatures.moment(x = mask,ref = DAPI)
  Cdh1_express <- computeFeatures.basic(x = mask,ref = Cdh1)
  Tubb3_express <- computeFeatures.basic(x = mask,ref = Tubb3)
  Zeb1_express <- computeFeatures.basic(x = mask,ref = Zeb1)
  
  tmp <- Result[[i]]
  rownames(tmp) <- tmp$cell
  tmp$x <- loc[,1]
  tmp$y <- loc[,2]
  tmp$Cdh1 <- Cdh1_express[,1]
  tmp$Tubb3 <- Tubb3_express[,1]
  tmp$Zeb1 <- Zeb1_express[,1]
  
  # generate the nuclear region
  whitePixelMask <- thresh(filter2(DAPI,gaussian_kernel(3,15)),w=15, h=15, offset=0.03)
  brush <- makeBrush(5, shape = "diamond")
  brushed_mask<- opening(whitePixelMask,kern = brush) %>% fillHull()
  
  nmask <- brushed_mask*mask # generate the nuclear mask

  Zeb1_nuclear <- computeFeatures.basic(nmask,Zeb1) %>% as.data.frame()
  Zeb1_nuclear$cell <- rownames(Zeb1_nuclear)
  
  tmp$Zeb1_nuclear <- NA
  tmp[Zeb1_nuclear$cell,"Zeb1_nuclear"] <- Zeb1_nuclear$b.mean
  
  Result[[i]] <- tmp
  
  rm(DAPI,Cdh1,Tubb3,Zeb1,mask,nmask,Cdh1_express,Tubb3_express,Zeb1_express,Zeb1_nuclear,tmp,loc,sample)
}
Origin <- data.table::rbindlist(Result) %>% as.data.frame()
Tumor <- Origin[which(Origin$annotation=="tumor cells"),]

# Define cell subtypes
Tumor$subtype <- "Cdh1-/Tubb3-"
Tumor$subtype[which(Tumor$Cdh1>0.25 & Tumor$Tubb3<0.25)] <- "Cdh1+/Tubb3-"
Tumor$subtype[which(Tumor$Cdh1<0.25 & Tumor$Tubb3>0.25)] <- "Cdh1-/Tubb3+"
Tumor$subtype[which(Tumor$Cdh1>0.25 & Tumor$Tubb3>0.25)] <- "Cdh1+/Tubb3+"
subtype <- c("Cdh1-/Tubb3-"="#8491B4","Cdh1+/Tubb3-"="#FF6A6A","Cdh1-/Tubb3+"="#00A087","Cdh1+/Tubb3+"="#B09C85")
Tumor$mouse <- str_extract(Tumor$sample,"C\\d+")
Tumor$Site <- str_extract(Tumor$sample,"(C\\d+-)(.*)(_Series)",group = 2)
save(Result,Origin,Tumor,subtype,file = "Result/Metadata.RData")


load("Metadata.RData")
plotdata <- Tumor[which(Tumor$subtype %in% c("Cdh1+/Tubb3-","Cdh1-/Tubb3+")),]

pdf("Nuclear_Zeb1_expression.pdf",height = 4,width = 4)
ggplot(plotdata,aes(subtype,Zeb1))+
  geom_violin(aes(fill=`subtype`),scale = "width",show.legend = F)+scale_x_discrete(limits=rev)+
  geom_boxplot(width=0.1,outliers = F)+scale_fill_manual(values = subtype)+theme_classic()+
  stat_compare_means(comparisons = list(c("Cdh1+/Tubb3-","Cdh1-/Tubb3+")),method = "wilcox.test")+coord_cartesian(clip = "off")
dev.off()
