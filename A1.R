setwd("~/Desktop/Spatial Ecology/Week_1")

#Install packages
install.packages(c("terra","sf","spatstat","cowplot","ggplot2","dismo","glmnet","maxnet","precrec","mlr","mapview","randomForest","rnaturalearth","dplyr"))
library(terra) #for spatial data
library(sf) #data frames with geometry
library(spatstat) #for point process modelling and converting between raster and pixel image objects
library(dismo)

# loading in my raster data
meles<- read.csv("Melesmeles.csv")



