setwd("~/Desktop/Spatial Ecology")

#Install packages
install.packages(c("terra","sf","mapview"))
library(terra)
library(sf)

# loading in my raster data
meles<- read.csv("Melesmeles.csv")