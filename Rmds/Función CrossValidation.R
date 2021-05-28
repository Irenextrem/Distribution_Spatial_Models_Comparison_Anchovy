################################################################################
#### CROSS VALIDATION PARA CAPACIDAD PREDICTIVA DE MODELOS EN INLA ####
################################################################################

#Librerias
library(sp)
library(rgdal) 
library(carData)
library(car)
library(nlme)
library(gstat)
library(sf)
library(spData)
library(spdep)
library(lattice)
library(survival)
library(Formula)
library(Hmisc)
library(raster) 
library(leaflet)
library(GGally)
library(maptools)
library(rgeos)
library(maptools) 
library(parallel)
library(foreach)
library(INLA) 
library(dotCall64)
library(grid)
library(spam)
library(fields)

library(dplyr)
library(PresenceAbsence)
library(MatrixModels)


#anchoa 
sdmdata <- read.csv("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Datos/sdmdata_atln_modelos.csv")
sdmdata <- sdmdata[,c(-1)]


#shapefile contiene los mapas unidos  de España, Francia y Portugal
paises <- readRDS("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Rmd/mapa_recortado.rds") #Mapa final recortado


#recortar zona de interés 
ext<-extent(-11,9,43,60)
cat <- crop(paises, ext) 
formula0 <- y~ -1 + beta0  + bathy + odismean +tempmean + salinity
formula1 <- y ~ -1 + beta0  + bathy + odismean +tempmean + salinity+  f(spatial,model=spde)

#cargar predictores
files<-(list.files("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Predictores/atl/Seleccionadas", full.names=T, pattern=".tif"))
predictors <- stack(files)
names(predictors) <- c("bathy","odismean","salinity","tempmean")


###########################################################################
### FUNCION DE VALIDACION CRUZADA 10 VECES ###
###########################################################################

simulacro <- function(DATAANCHOA, CAT, FORMULA, PREDICTORS){
  
  #### Dividir mi banco de datos en taining y test ####
  
  splitdf <- function(dataframe,fraction, seed=NULL) {
    
    if (!is.null(seed)) set.seed(seed)
    
    index <- 1:nrow(dataframe)
    
    trainindex <- sample(index, trunc(length(index)*fraction))
    
    trainset <- dataframe[trainindex, ]
    
    testset <- dataframe[-trainindex, ]
    
    list(trainset=trainset,testset=testset)
    
  }
  
  splits <- splitdf( DATAANCHOA,    #DATAANCHOA,
                     0.8, seed= i )
  
  training <- splits$trainset
  test <- splits$testset
  
  
  ##### Ajustar el modelo espacial  con el trainset #####
  
  ### -- Definir polygono y crear MESH -- ###
  
  xym <- as.matrix(data.frame(x=c(-11,0,0,9,9,
                                  -11,-11),
                              y=c(43,43,43,55,60,
                                  60,43)))
  
  p = Polygon(xym)
  ps = Polygons(list(p),1) 
  sps = SpatialPolygons(list(ps))
  
  map_rec<-crop( CAT, #CAT,
                 sps)
  proj4string(sps)<-proj4string(map_rec)
  coast <- gDifference(sps, map_rec )
  
  ### -- Definir MESH -- ###
  boundary<-inla.nonconvex.hull(as.matrix(training[,6:7]))
  mesh<-inla.mesh.2d(boundary=boundary, max.edge=c(1, 3), 
                     cutoff=0.01,  offset=c(-0.9, -0.09))
  
  ### -- Definir SPDE -- ###
  spde <- inla.spde2.matern(mesh)
  
  ### -- Matriz que une los datos con el MESH -- ###
  A.est <- inla.spde.make.A(mesh, loc=cbind(training$x, training$y))
  
  ### -- inla.stack para estimar -- ###
  Stk.est<-inla.stack(data=list(y=training$pb), #$ anchoa ; $sardina (cambiar dependiendo del banco usado)
                      A=list(A.est, 1),
                      effects=list(spatial=1:spde$n.spde,
                                   data.frame(beta0=1, training[,2:5])),
                      tag='est')
  
  ### -- Estimar el modelo con training -- ###
  
  model.est <- inla(FORMULA, #Poner formula0 y formula 
                    data=inla.stack.data(Stk.est), family="binomial" ,
                    control.compute=list(dic=TRUE,cpo=TRUE, waic=TRUE, return.marginals=TRUE), 
                    control.predictor=list(A=inla.stack.A(Stk.est), compute=TRUE, 
                                           quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975)),     
                    num.threads = 4,
                    verbose=F)
  
  
  
  ##### Hacer la predicción con el trainset #####
  
  
  ### -- Definir SPDE -- ###
  spde <- inla.spde2.matern(mesh)
  
  ### -- Definir que la predicción se haga SOLO dentro del mar -- ###
  dxy <- apply(bbox(coast),1, diff)
  r <- dxy[1]/dxy[2]
  m<-60
  proj.grid.mat <- inla.mesh.projector(mesh, 
                                       xlim=bbox(coast)[1,],
                                       ylim=bbox(coast)[2,] ,
                                       dims=c(r, 1)*m)
  
  ov <- over(SpatialPoints(proj.grid.mat$lattice$loc, coast@proj4string), 
             coast) #hago una limpieza definiendo NAs los valores que se encuentran fuera de mi boundary
  
  i.map <- is.na(ov) # chequeo los puntos del grid que esten dentro de mi mapa
  
  ### --- Matrix proyectará las coordenadas del test al mesh--- ###
  a.pred <- inla.spde.make.A(mesh, loc=proj.grid.mat$lattice$loc[!i.map, ])
  
  ### --- Stack de la predicción --- ###
  stk.pred <- inla.stack(data=list(y=NA),
                         A=list(a.pred, 1), 
                         effects=list(spatial=1:spde$n.spde,
                                      data.frame(beta0 = 1, 
                                                 extract( PREDICTORS, #PREDICTORS, 
                                                          proj.grid.mat$lattice$loc[!i.map, ]))),
                         tag='pred')
  
  stk <- inla.stack(Stk.est, stk.pred)
  
  
  ### --- Modelo Predictivo --- ###
  
  model.pred <- inla(FORMULA, 
                     data=inla.stack.data(stk), family="binomial",
                     control.predictor=list(A=inla.stack.A(stk), compute=TRUE, link=1), 
                     control.inla=list(strategy = "simplified.laplace"), # Strategy
                     #control.mode=list(model.est$mode$theta, restart=TRUE), #Mode 
                     control.results=list(return.marginals.random=FALSE,
                                          return.marginals.predictor=FALSE), # Avoid some marginals
                     num.threads = 4,
                     verbose=FALSE) # La predicción se hace con training
  
  #Error in inla.check.control(control.mode, data) : 
  #Name `' in control-argument `control.mode', is void.
  #  Valid ones are:
  #	fixed
  #	restart
  #	result
  #	theta
  #	x
  
  
  ### --- Indice para prediccion del test --- ###
  idx <- inla.stack.index(stk, 'pred')$data #extraer los valores de la prediccion con un indice
  
  
  
  ##### Crear dataset unico co presencias/ausencias para mis valores actuales y predichos #####
  
  # Hay que hacer un unico dataframe con coordenadas, presencia/ausencia del test y
  #Prob Media de presencia/ausencia (en 0/1) del modelo predictivo
  
  prob.mean  <- matrix(NA, proj.grid.mat$lattice$dims[1], proj.grid.mat$lattice$dims[2]) #Calculo la media
  prob.mean[!i.map] <- c(model.pred$summary.fitted.values$mean[idx])
  
  prob.mean.raster<-raster(list(x = proj.grid.mat$x, 
                                y = proj.grid.mat$y,
                                z = prob.mean ))
  
  coords=cbind(test$x,test$y) #Coordenadas de los datos del test
  bb<-extract(prob.mean.raster,coords) #Extraigo del raster los valores para las coordenadas
  
  #test$prob_mean<-as.data.frame(bb, xy=TRUE)
  
  prueba<-as.data.frame(cbind(coords,bb)) #Lo metemos en una matriz
  colnames(prueba)<-c("lon", "lat", "mean_predic")
  #prueba$predict<- as.integer(ifelse(prueba$mean_predic >= 0.7 , "1", "0"))
  
  m3=merge(test,prueba, by.x=c("x","y"),by.y=c("lon", "lat")) #, all=T //Uno valores de la prueba con el test
  
  ### -- Mi dataframe final -- ###
  ID <-  as.numeric(dimnames(m3)[[1]])  # Identificador
  verif<-as.data.frame(cbind(ID,m3$pb,m3$mean_predic))  #los uno en una tabla en forma: ID, Observados, predichos
  
  
  ##### MIRAR AUC y TSS #####
  
  
  Threshold<-optimal.thresholds(verif,threshold=101,
                                which.model=1:(ncol(verif)-2), #which.model=1:(ncol(verif)-2)
                                opt.methods=4, na.rm=T)
  
  confusionmatrix <-cmx(verif,threshold=Threshold[,2], which.model = 1, na.rm=T)
  
  #aucc<- auc(confusionmatrix, which.model = 1, na.rm=T) da error y no lo hace
  
  aucc<- auc(verif, which.model = 1, na.rm=T)
  #auc_mean<- aucc$AUC
  #☻auc_sd<- aucc$AUC.sd
  
  kappa<-Kappa(confusionmatrix)
  #kappa_mean<- kappa$Kappa
  #kappa_sd<- kappa$Kappa.sd
  
  
  sensi<- sensitivity(confusionmatrix)
  #sensi_mean<- sensi$sensitivity
  #sensi_sd<- sensi$sensitivity.sd
  
  speci<-specificity(confusionmatrix)
  #speci_mean<- speci$specificity
  #speci_sd<- speci$specificity.sd
  
  
  #### Me devuelve los valores calculados en un dataframe ####
  
  return(c(aucc$AUC, aucc$AUC.sd, kappa$Kappa, kappa$Kappa.sd, sensi$sensitivity, sensi$sensitivity.sd, speci$specificity, speci$specificity.sd))
  #return(aucc, kappa, sensi, speci)
  
  
}

# calcular el TSS 

cosinas_ire <- matrix(ncol=8,nrow=10)
cosinas_ire1 <- matrix(ncol=8,nrow=10)


for(i in 1:10){cosinas_ire[i,]<- simulacro(sdmdata, cat, formula0, predictors3)} 
for(i in 1:10){cosinas_ire1[i,]<- simulacro(sdmdata, cat, formula1, predictors3)} 

cosinas_ire3 <- cbind(aucc$AUC, aucc$AUC.sd, kappa$Kappa, kappa$Kappa.sd, sensi$sensitivity, sensi$sensitivity.sd, speci$specificity, speci$specificity.sd)
cosinas_ire31 <- cbind(aucc$AUC, aucc$AUC.sd, kappa$Kappa, kappa$Kappa.sd, sensi$sensitivity, sensi$sensitivity.sd, speci$specificity, speci$specificity.sd)

cosinas0<-cosinas_ire 
cosinas1<-cosinas_ire1

colnames(cosinas0) <- c('AUC_mean','AUC_sd','Kappa_mean','Kappa_sd', 'Sensi_mean','Sensi_sd','Speci_mean','Speci_sd')
colnames(cosinas1) <- c('AUC_mean','AUC_sd','Kappa_mean','Kappa_sd', 'Sensi_mean','Sensi_sd','Speci_mean','Speci_sd')

cosinas0<-apply(cosinas0, 2, mean)
cosinas1<-apply(cosinas1, 2, mean)

write.csv(cosinas0,"C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/GLM INLA/cv_inla_sin_es.csv")
write.csv(cosinas1,"C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/GLM INLA/cv_inla_es.csv")

CVanchoa<-read.csv("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/GLM INLA/cv_inla_sin_es.csv")
CVanchoa<-CVanchoa[,-1]
CVanchoa<- matrix(CVanchoa,nrow = 1)

CVanchoa1<-read.csv("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/GLM INLA/cv_inla_es.csv")
CVanchoa1<-CVanchoa1[,-1]
CVanchoa1<- matrix(CVanchoa1,nrow = 1)

colnames(CVanchoa) <- c('AUC_mean','AUC_sd','Kappa_mean','Kappa_sd', 'Sensi_mean','Sensi_sd','Speci_mean','Speci_sd')
colnames(CVanchoa1) <- c('AUC_mean','AUC_sd','Kappa_mean','Kappa_sd', 'Sensi_mean','Sensi_sd','Speci_mean','Speci_sd')

TSS0 <-CVanchoa[,5] + CVanchoa[,7] -1 #Hago nueva columna para calcular el TSS
CV0 <- matrix(cbind(TSS0,CVanchoa),nrow = 1,ncol = 9)
colnames(CV0) <- c("TSS",'AUC_mean','AUC_sd','Kappa_mean','Kappa_sd', 'Sensi_mean','Sensi_sd','Speci_mean','Speci_sd')


TSS1<-CVanchoa1[,5] + CVanchoa1[,7] -1 #Hago nueva columna para calcular el TSS
CV1 <- matrix(cbind(TSS1,CVanchoa1),nrow = 1,ncol = 9)
colnames(CV1) <- c("TSS",'AUC_mean','AUC_sd','Kappa_mean','Kappa_sd', 'Sensi_mean','Sensi_sd','Speci_mean','Speci_sd')

write.csv(CV0,"C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/GLM INLA/cv_inla_sin_estss.csv")
write.csv(CV1,"C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/GLM INLA/cv_inla_estss.csv")
