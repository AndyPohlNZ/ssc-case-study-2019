rm(list = ls())
library(raster)
library(rgdal)
library(plotly)
library(Thermimage)

### label file ####
setwd("~/Downloads/SSC19-master")#please change the directory to your working path

train_label <- read.csv("train_label.csv",header = T,row.names = 1)
test_label <- read.csv("test_label.csv",header = T,row.names = 1)

### read training images ####
setwd("~/Downloads/SSC19-master/train")#please change the directory to your working path
##create an array (train_array) containing pixel value of 2400 training images
train_array<-array(NA,c(520,696,2400))
#start_time <- Sys.time()
for (i in 1:2400){
  train_array[,,i]<-as.array(raster(as.character(train_label$image_name[i])))
}#Large array, 3.2 Gb
#end_time <- Sys.time()
#end_time - start_time #Time difference of 10.15881 mins
##plot 1 image as an example
#plot_ly(x=NULL,y=NULL,z=flip.matrix(train_array[,,1]), type = "heatmap",colors=colorRamp(c("black", "white")) ,showscale = F)

### read testing images ####
setwd("~/Downloads/SSC19-master/test")#please change the directory to your working path
##create an array (test_array) containing pixel value of 1200 testing images
test_array<-array(NA,c(520,696,1200))
#start_time <- Sys.time()
for (i in 1:1200){
  test_array[,,i]<-as.array(raster(as.character(test_label$image_name[i])))
}
#end_time <- Sys.time()
#end_time - start_time #Time difference of 5.576581 mins



