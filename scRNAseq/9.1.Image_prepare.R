library(jpeg)
library(EBImage)
library(stringr)
library(dplyr)

# Input path with each IF image separated in single channel, appendixed as _ch01, _ch02, etc.:
file_list <- list.files("Image/")
dir.create("Merge")

for (Sample in file_list) {
  roi <- list.files(paste0("Image/",Sample))
  for (r in roi) {
    print(paste(Sample,r,sep = "_"))
    DAPI <- readImage(paste0("Image/",Sample,"/",r,"/",r,"_ch03.tif"))
    Tubb3 <- readImage(paste0("Image/",Sample,"/",r,"/",r,"_ch02.tif"))
    Cdh1 <- readImage(paste0("Image/",Sample,"/",r,"/",r,"_ch01.tif"))
    img <- rgbImage(red = Cdh1,green = Tubb3, blue = DAPI)
    writeImage(img,files = paste0("Merge/",Sample,"_",r,".tif")) # Merge images without Zeb1 channel and output to Merge/
    rm(DAPI,Cdh1,Tubb3,img)
  }
}

