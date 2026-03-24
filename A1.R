setwd("~/Desktop/Spatial Ecology/Week_1+3data")

#----------------------------installing packages--------------------------------

install.packages(c("terra","sf","spatstat","cowplot","ggplot2","dismo","glmnet","maxnet","precrec","mlr","ranger","mapview"))
library(terra) #for spatial data
library(sf) #data frames with geometry
library(spatstat) #for point process modelling and converting between raster and pixel image objects
library(dismo) 
library(mlr)
library(ranger)


#------------------------loading and cleaning my data---------------------------

#Function that will later convert raster objects to pixel images
raster.as.im = function(im) {
  r = raster::res(im)[1]
  orig = ext(im)[c(1,3)]
  xx = orig[1] + seq(from=0,to=(ncol(im) - 1)*100,by=r)
  yy = orig[2] + seq(from=0,to=(nrow(im) - 1)*100,by=r)
  
  #building a single matrix
  mat=matrix(raster::values(im), ncol = ncol(im), 
             nrow = nrow(im), byrow = TRUE)[nrow(im):1, ]
  return(spatstat.geom::im(mat, 
                           xcol = xx, yrow = yy))
}

#loading in my data
meles<- read.csv("Melesmeles.csv")
#subset to cases with coordinate uncertainty <100 metres
meles<-meles[meles$Coordinate.uncertainty_m<=1000,]

#making spatial points layer and crs object
meles.latlong=data.frame(x=meles$Longitude,y=meles$Latitude)
meles.sp=st_as_sf(meles.latlong,coords=c("x","y"),crs="epsg:4326")

#loading in study area (Scotland) and cropping
scot=st_read('scotSamp.shp')
LCM=rast("LCMUK.tif")
LCM=crop(LCM,st_buffer(scot, dist= 1000))

#aggregate data and increase cell size
LCM=aggregate(LCM$LCMUK_1,fact=4,fun="modal")

#Making sure that meles points are the same projection as the LCM layer and crop
meles.sp=st_transform(meles.sp,crs(LCM))
melesFin=meles.sp[scot,]
LCM=crop(LCM,scot,mask=TRUE)

#inspect
plot(LCM)
plot(melesFin$geometry,add=T)

#reclassify the land cover raster (covariates)
LCM=as.factor(LCM$LCMUK_1)
reclass = c(0,1,rep(0,20))
RCmatrix=cbind(levels(LCM)[[1]],reclass)
RCmatrix=RCmatrix[,2:3]

#making sure new columns are numeric 
RCmatrix=apply(RCmatrix,2,FUN=as.numeric)
broadleaf=classify(LCM, RCmatrix)

#generating a raster layer for the proportion of broadleaf woodland within a 1800 metre radius.
nPix=round(1800/res(LCM)[1])
nPix=(nPix*2)+1
weightsMatrix=matrix(1:nPix^2,nrow=nPix,ncol=nPix)

#get focal cell 
x=ceiling(ncol(weightsMatrix)/2)
y=ceiling(nrow(weightsMatrix)/2)
focalCell=weightsMatrix[x,y]
indFocal=which(weightsMatrix==focalCell,arr.ind = TRUE)

#populating matrix with distances (limit is 1800 metres)
distances=list()
for(i in 1:nPix^2){
  ind.i=which(weightsMatrix==i,arr.ind=T)
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(LCM)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(LCM)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distances[[i]]=dist.i
  
}

#entering values into the matrix
weightsMatrix[]=unlist(distances)
weightsMatrix[weightsMatrix>1800]=NA

#plot weights matrix
plot(rast(weightsMatrix))

#normalising the weights matrix
weightsMatrixNorm=weightsMatrix
weightsMatrixNorm[!is.na(weightsMatrixNorm)]=1/length(weightsMatrixNorm[!is.na(weightsMatrixNorm)])
sum(weightsMatrixNorm,na.rm=T)
plot(rast(weightsMatrixNorm))

#apply weights matrix to the woodland layer
lcm_wood_1800=focal(broadleaf,w=weightsMatrixNorm,fun="sum")
plot(lcm_wood_1800)

#creating layer for urban class
reclassUrban = c(rep(0,19),1,1)
RCmatrixUrban= cbind(levels(LCM)[[1]],reclass)
RCmatrixUrban=RCmatrixUrban[,2:3]

#making sure new columns are numeric
RCmatrixUrban=apply(RCmatrixUrban,2,FUN=as.numeric)
urban = classify(LCM, RCmatrixUrban)

#get number of pixcels needed to cover the 2300 metre radius for the urban class
nPixUrban=round(2300/res(LCM)[1])
nPixUrban=(nPixUrban*2)+1

#building weights matrix
weightsMatrixUrban=matrix(1:nPixUrban^2,nrow=nPixUrban,ncol=nPixUrban)

#get focal cell 
x=ceiling(ncol(weightsMatrixUrban)/2)
y=ceiling(nrow(weightsMatrixUrban)/2)
focalCell=weightsMatrixUrban[x,y]
indFocal=which(weightsMatrixUrban==focalCell,arr.ind = TRUE)

#computing distances
distancesUrban=list()
for(i in 1:nPixUrban^2){
  ind.i=which(weightsMatrixUrban==i,arr.ind=T)
  diffX=abs(ind.i[1,1]-indFocal[1,1])*res(LCM)[1]
  diffY=abs(ind.i[1,2]-indFocal[1,2])*res(LCM)[1]
  
  dist.i=sqrt(diffX^2+diffY^2)
  distancesUrban[[i]]=dist.i
  
}

#weights matrix
weightsMatrixUrban[]=unlist(distancesUrban)
weightsMatrixUrban[weightsMatrixUrban>2300]=NA
weightsMatrixUrban[!is.na(weightsMatrixUrban)]=1/length(weightsMatrixUrban[!is.na(weightsMatrixUrban)])

#sum urban class from all surrounding cells
lcm_urban_2300=focal(urban,w=weightsMatrixUrban,fun="sum")

#considering elevation
demScot=rast('demScotland.tif')
demScot=terra::resample(demScot,lcm_wood_1800)
plot(demScot)

#stacking the covariate layers together
allEnv=c(lcm_wood_1800,lcm_urban_2300,demScot)
names(allEnv)=c("broadleaf","urban","elev")

#Creating background points
set.seed(11)
back = spatSample(allEnv,size=2000,as.points=TRUE,method="random",na.rm=TRUE) 
back=back[!is.na(back$broadleaf),]
back=st_as_sf(back,crs="EPSG:27700")

#extraction of data to points
eP=terra::extract(allEnv,melesFin)

#binding presence data 
Pres.cov=st_as_sf(cbind(eP,melesFin))
Pres.cov$Pres=1

#Removing first column
Pres.cov=Pres.cov[,-1]

#getting and combining coordinates
coordsPres=st_coordinates(Pres.cov)
Back.cov=st_as_sf(data.frame(back,Pres=0))
coordsBack=st_coordinates(back)
coords=data.frame(rbind(coordsPres,coordsBack))
colnames(coords)=c("x","y")

#mlr
all.cov=rbind(Pres.cov,Back.cov)
all.cov=cbind(all.cov,coords)
all.cov=na.omit(all.cov)
all.cov=st_drop_geometry(all.cov)

#making target variable categorical
task=all.cov
head(all.cov)
task$Pres=as.factor(task$Pres)
task = makeClassifTask(data = task[,c(1:4)], target = "Pres",
                       positive = "1", coordinates = task[,5:6])



#-----------------------Binomial--------------------------
lrnBinomial = makeLearner("classif.binomial",
                          predict.type = "prob",
                          fix.factors.prediction = TRUE)


#five folds 
perf_levelCV = makeResampleDesc(method = "RepCV", predict = "test", folds = 5, reps = 5)

#spatial cross-validation
perf_level_spCV = makeResampleDesc(method = "SpRepCV", folds = 5, reps = 5)
cvBinomial = resample(learner = lrnBinomial, task =task,
                      resampling = perf_levelCV, 
                      measures = mlr::auc,
                      show.info = FALSE)

print(cvBinomial)


#plot locations
plots = createSpatialResamplingPlots(task,resample=cvBinomial,
                                     crs=crs(allEnv),datum=crs(allEnv),color.test = "red",point.size = 1)
library(cowplot)
cowplot::plot_grid(plotlist = plots[["Plots"]], ncol = 3, nrow = 2,
                   labels = plots[["Labels"]])

#GRAPH 1: Random 5-FOLD 

#Carrying out and plotting Binomial spatial cross validation
sp_cvBinomial = resample(learner = lrnBinomial, task =task,
                         resampling = perf_level_spCV, 
                         measures = mlr::auc,
                         show.info = FALSE)

print(sp_cvBinomial)

#AUC value: 0.7598611

plotsSP = createSpatialResamplingPlots(task,resample=sp_cvBinomial,
                                       crs=crs(allEnv),datum=crs(allEnv),color.test = "red",point.size = 1)
cowplot::plot_grid(plotlist = plotsSP[["Plots"]], ncol = 3, nrow = 2,
                   labels = plotsSP[["Labels"]])

#GRAPH 2: Spatial 5-fold repeated CV


#------------------Relationships between covariates and Meles-----------------------

#binomal glm.
glm.meles=glm(Pres~broadleaf+urban+elev,binomial(link='logit'),
                data=all.cov)
prGLM=predict(allEnv,glm.meles,type="response")
plot(prGLM)

#MAP 1:COMBINED COVARIATES AND MELES PREDICTION OCCURANCE

#INDEPENDANT COVARIANT: broadleaf
#building new data frame 
glmNew=data.frame(broadleaf=seq(0,max(all.cov$broadleaf),length=1000),
                  elev=mean(all.cov$elev),
                  urban=mean(all.cov$urban))

preds = predict(glm.meles, newdata = glmNew, type = "response", se.fit = TRUE)
glmNew$fit = preds$fit
glmNew$se = preds$se.fit

head(glmNew)

#Plotting broadleaf
library(ggplot2)
ggplot(glmNew, aes(x = broadleaf, y = fit)) +
  
  geom_ribbon(data = glmNew, aes(y = fit, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "green", alpha = 0.3) +
  geom_line(data = glmNew, aes(y = fit)) 

#GRAPH 3:Broadleaf and badger response curve


#INDEPENDANT COVARIANT: Urban
#building new data frame 
glmNewUrban=data.frame(urban=seq(0,max(all.cov$urban),length=1000),
                       elev=mean(all.cov$elev),
                       broadleaf=mean(all.cov$broadleaf))

predUrban = predict(glm.meles, newdata = glmNewUrban, type = "response", se.fit = TRUE)
glmNewUrban$fit = predUrban$fit
glmNewUrban$se = predUrban$se.fit

#Plotting Urban
ggplot(glmNewUrban, aes(x = urban, y = fit)) +
  
  geom_ribbon(data = glmNewUrban, aes(y = fit, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "red", alpha = 0.3) +
  geom_line(data = glmNewUrban, aes(y = fit)) 

#GRAPH 4: Urban and badger response curve

#INDEPENDANT COVARIANT:elevation
#building new data frame 
glmNewElev=data.frame(elev=seq(0,max(all.cov$elev),length=1000),
                      urban=mean(all.cov$urban),
                      broadleaf=mean(all.cov$broadleaf))

predElev = predict(glm.meles, newdata = glmNewElev, type = "response", se.fit = TRUE)
glmNewElev$fit = predElev$fit
glmNewElev$se = predElev$se.fit

#Plotting elevation
ggplot(glmNewElev, aes(x = elev, y = fit)) +
  
  geom_ribbon(data = glmNewElev, aes(y = fit, ymin = fit - 1.96 * se, ymax = fit + 1.96 * se),
              fill = "blue", alpha = 0.3) +
  geom_line(data = glmNewElev, aes(y = fit)) 

#GRAPH 5: Elevation and badger response curve





#----------------------------Point Process Model---------------------------------
#converting our raster covariates to image objects.
broadleafIm=raster.as.im(raster(allEnv$broadleaf))
urbanIm=raster.as.im(raster(allEnv$urban))
elevIm=raster.as.im(raster(allEnv$elev))
window.poly=as.owin(urbanIm)

#Setting of analysis window
melesCoords=st_coordinates(melesFin) 
pppmeles=ppp(melesCoords[,1],melesCoords[,2], window = window.poly)

#plot
plot(allEnv$broadleaf)
plot(pppmeles,add=T)
pppmeles=as.ppp(pppmeles)

#Rescaling data from m to km
pppmeles=rescale(pppmeles, 1000)
broadleafIm=rescale(broadleafIm,1000)
elevIm=rescale(elevIm,1000)
urbanIm=rescale(urbanIm,1000)

#Ripley’s K test
Kcsr=envelope(pppmeles,Kest,nsim=39,VARIANCE=T,nSD=1,global =TRUE)
plot(Kcsr,shade=c("hi","lo"),legend=T)

#GRAPH 6: Ripley's K-function

#Using AIC to create a better model 
ndTry=seq(100,1000,by=100)
for(i in ndTry){
  Q.i=quadscheme(pppmeles,method = "grid",nd=i)
  fit.i=ppm(Q.i~broadleafIm+elevIm+urbanIm)
  print(i)
  print(AIC(fit.i))
}
Q=quadscheme(pppmeles,method = "grid",nd=900)

#plot for each covariate
plot(rhohat(pppmeles,broadleafIm)) #broadleaf --------GRAPH 7
plot(rhohat(pppmeles,urbanIm)) #urban -------------GRAPH 8
plot(rhohat(pppmeles,elevIm)) #elevation --------------GRAPH 9


#THRID ORDER POLYNOMIAL TREND
firstPPMod=ppm(Q~poly(broadleafIm,3)+poly(elevIm,2)+poly(urbanIm,2)+x+y)

#---------------------------diagnostic tests-----------------------

#influenced by spatial dependence.
firstModEnv=envelope(firstPPMod,Kest,nsim=39,VARIANCE=TRUE,nSD=1,global =TRUE)
plot(firstModEnv)

#GRAPH 10: Spatial randomness

#Thomas process
thomasMod=kppm(Q~poly(broadleafIm,3)+poly(elevIm,2)+poly(urbanIm,2)+x+y,"Thomas")

#simulate inhomogeneous patterns).
thomasEnv=envelope(thomasMod,Kinhom,nsim=39,VARIANCE=TRUE,nSD=1,global=TRUE)
plot(thomasEnv)

#GRAPH 11: Thomas process


#AUC
plot(roc(thomasMod))

#Graph 12: AUC value for PPM

#Print
auc.kppm(thomasMod)

#mapped output
prPPMod=predict(thomasMod)

plot(prPPMod)

#Map 2: Mapped results of PPM


