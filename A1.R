setwd("~/Desktop/Spatial Ecology/Week_1")

#Install packages
install.packages(c("terra","sf","mapview","dismo","glmnet","maxnet","randomForest","rnaturalearth","dplyr","precrec"))
library(terra)
library(sf)
library(dismo)

# loading in my raster data
meles<- read.csv("Melesmeles.csv")



