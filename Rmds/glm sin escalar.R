# Directorio de trabajo

setwd("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Rmd")

# Mapa del mundo
data(wrld_simpl) 

#Presencias
data <- read.csv('C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Datos/Anchoas_Aqua_atln.csv',TRUE,",")
data <- data[,-1] #Le quito al primera que es una columna de índices
colnames(data) <- c('Lon','Lat') #Le doy nombre a las columnas

#Predictores modelo
cuatro_pred <- (list.files("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Predictores/atl/Seleccionadas", full.names=T, pattern=".tif")) #Cargo bathy que lo he metido en una carpeta aparte
predictors3 <- stack(cuatro_pred)
names(predictors3) <- c("bathy","odismean","salinity","tempmean")

# Ausencias (Con todos los predictores pues es así como se han creado para meterlas al sdmdata_modelos)
files<-(list.files("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Predictores/atl", full.names=T, pattern=".tif"))#change directory
predictores <- stack(files)
names(predictores) <- c("bathy","chlomean","ppmean","odismean","salinity","tempmean")
predictors2 <- scale(predictores) #Los escalo
set.seed(141592) 
backgr <- randomPoints(predictors2, 1000) 
plot(predictors3)

# Base de datos sin escalar
sdmdata_nos <- read.csv("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Datos/sdmdata_sin_escalar.csv") #Sin las variables escaladas
sdmdata_nos<- sdmdata_nos[,c(-1)]

#GLM sin componente espacial
glm1 <- glm(pb ~ bathy + odismean + tempmean + salinity, data=sdmdata_nos, family="binomial") 

predm1 <- predict(predictors3,glm1,type='response')

# saveRDS(predm1, "C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/GLM frecuentista/Predicciones_glm_Atl.ascii")

ggplot() +
  geom_raster(data = raster::as.data.frame(predm1 , xy = TRUE) , aes(x = x, y = y, fill = layer)) +
  coord_equal() +
  labs( x = "", y = "")+theme_minimal()+
  scale_color_brewer(palette = 'YlGnBu')+ labs(fill='')


########    RAC   ###################
#Extract residuals from the glm model and map them
# r<- raster(xmn=-11, xmx=9, ymn=43, ymx=60, nrows=324, ncols=660);r #Hago un raster grande
# res(r) <‐ 0.083 #Le doy una resolución de 0.083
# xy <-cbind(sdmdata_nos$x, sdmdata_nos$y)#Cojo las coordenadas
# xy_residuals <-cbind(xy, resid(glm1)) #Uno las coordenadas y los residuos a un mismo objeto
# par(mfrow=c(1,2))
# r[cellFromXY(r, xy_residuals)] <-xy_residuals[,3] #Doy esos valoes de residuos a esas coordenadas
# plot(r,col='red') #Efectivamente me salen residuos
# ext <- c(-11,9,43,60) #Extensión a cortar
# salinity <- raster("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Predictores/RAC/GLM/salinity_atln.TIF") #Cargo este raster para poder ajustar el raster r
# salinity
# 
# r <- resample(r,salinity) #Ajusto las dimensiones de r a las de salinity
# r<-crop(r,ext);r #Lo corto
# # writeRaster(r, filename="RAC_glm_atln.tif", format="GTiff", overwrite=TRUE) #Lo guado en el wd
# 
# 
# #Calculate residuals autocovariate
# focal_rac_rast <-focal(r, w=matrix(1,3,3), fun = mean,  na.rm = TRUE)
# 
# #Extract the values of the focal operation from focal_rac_rest raserfile using the coordinates stored in xy
# focal_rac_vect <-extract(focal_rac_rast, xy)
# 
# 
# #Add as a column to the data
# dd<-cbind(sdmdata_nos, focal_rac_vect)
# I <- is.na(dd$focal_rac_vect)
# ddd<- dd[!I,]
# dim(dd)
# dim(ddd)
# write.csv(ddd,'C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Datos/rac_glm_SINESC.csv')
ddd<- read.csv('C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Datos/rac_glm_SINESC.csv')

##############################
### --- Kriging con RAC--- ###
##############################

# salinity<- load_layers(c("BO_salinity"))# datadir = tempdir())
# e<- extent(-11,9,43,60)
# sal=crop(salinity,e) #Corto el raster
# 
# #Saco las coordenadas de uno de los rasters. Tengo que agregar porque sino no funciona el kriging.
# loci1 <-as.data.frame(coordinates(sal))
# summary(loci1)
# 
# x.range <- as.numeric(c(-10.958, 8.958))  # min/max longitude of the interpolation area
# y.range <- as.numeric(c(42.04,59.96))  #in/max latitude of the interpolation area
# grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by =  0.08333333), y = seq(from = y.range[1], 
#                                                                                           to = y.range[2], by =  0.08333333))  # expand points to grid
# coordinates(grd) <- ~x + y #Coordenadas del grid
# proj4string(grd) <- CRS("+init=epsg:4326") #Coordenadas geográficas
# gridded(grd)     <- TRUE  # Create SpatialPixel object
# fullgrid(grd)    <- TRUE  # Create SpatialGrid object
# 
# coords <- cbind(ddd$x,ddd$y) #Asigno coordenadas
# matrix<-cbind(ddd$focal_rac_vect,ddd[,c(6,7)]) #Uno las coordenadas a los valores de rac
# matrix<-as.data.frame(matrix) #Lo paso a data.frame
# colnames(matrix)<-c("Pred","Lon","Lat") #Renombro
# 
# coordinates(matrix) <- c("Lon", "Lat")  #Doy nombre a las coordenadas
# proj4string(matrix) <- CRS("+init=epsg:4326") #Coordenadas geográficas
# 
# #Interpolate the grid cells using a power value of 2 (idp=2.0)
# rac <- gstat::idw(Pred~ 1, locations = matrix, newdata=grd, idp=0.1) 
# 
# # Plot
# plot(rac,col=tim.colors(100)[1:100],main=" ", axes=T)
# data(wrld_simpl)
# plot(wrld_simpl, add=TRUE,col='dark grey')
# 
# #convert to raster
# RAC=as.data.frame(rac) #Paso la interpolación a data.frame
# coordinates(RAC) <- ~ x + y #Le asigno unas coordenadas
# gridded(RAC) <- TRUE #Creo un objeto pixelado espacial
# RAC <- raster(RAC) #Lo transformo en raster
# 
# plot(RAC)
# 
# #####################################################
# #  CUT WITH BATHY
# ##################################################
# 
# depth<- load_layers(c("MS_bathy_5m")) #Cargo batimetría
# depth=crop(depth,e) #La corto
# depth=abs(depth) ##funcao abs converte em positivo ### formato raster
# 
# ## tem que transformar em data.frame
# 
# matrix<- cbind(coordinates(depth), depth=getValues(depth)) #Cojo los valores de batimetría para esas coordenadas
# I <- is.na(matrix[,3]) #Identifico NAs
# matrix<- matrix[!I,] #Quito los NAs
# matrix<-as.data.frame(matrix) #Transformo a data.frame
# new<- subset(matrix, matrix$depth > 0) ###trocar profundidade
# 
# ###transformar data.frame em raster de novo
# 
# xy <- cbind(new$x, new$y) #Me quedo con las coordenadas
# rast<- raster(xmn=-20, xmx=35, ymn=43, ymx=70, nrows=276, ncols=624)#2  
# p<- rasterize(xy, rast, new$depth, fun=max,na.rm=F) #Hago un raster
# p<-resample(p,sal) #Le doy las dimensiones de salinidad
# e<-extent(p) #
# 
# sp<-crop(RAC,e)
# sp=resample(sp,p)
# sp<-raster::mask(sp,p)
# 
# 
# # Plot
# plot(sp,col=tim.colors(100)[1:100],main=" ", axes=T)
# data(wrld_simpl)
# plot(wrld_simpl, add=TRUE,col='dark grey')
# 
# plot(RAC,col=tim.colors(100)[1:100],main=" ", axes=T)
# data(wrld_simpl)
# plot(wrld_simpl, add=TRUE,col='dark grey')
# 
# writeRaster(sp, filename="sp_glm_atln_sinescalar.tif", format="GTiff", overwrite=TRUE)

# Predictores con RAC
files<-(list.files("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Predictores/RAC/GLM", full.names=T, pattern=".tif"))
predictors <- stack(files)
names(predictors) <- c("bathy","odismean","salinity","focal_rac_vect","tempmean")
predictors_rac <- scale(predictors) 

#Modelo con efecto espacial
glm_rac <-  glm(pb ~ bathy + odismean + tempmean + salinity + focal_rac_vect, data=ddd, family="binomial") 

# PREDICCIÓN
predm_glm_rac<-predict(predictors_rac, glm_rac,type='response')

saveRDS(predm_glm_rac, "C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/GLM frecuentista/Predicciones_glm_rac_Atl_sinescalar.ascii")

ggplot() +
  geom_raster(data = raster::as.data.frame(predm_glm_rac , xy = TRUE) , aes(x = x, y = y, fill = layer)) +
  coord_equal() +
  labs( x = "", y = "")+theme_minimal()+
  scale_color_brewer(palette = 'YlGnBu')+ labs(fill='')

plot(predm_glm_rac, col=tim.colors(100)[1:100],main="Glm RAC", axes=T)
data(wrld_simpl)
plot(wrld_simpl, add=TRUE,col='dark grey')

# Plots 

par(mfrow=c(1,2))

## PREDICTIVA SIN RAC
plot(predm1, col=tim.colors(100)[1:100],main="Glm sin RAC", axes=T,zlim=c(0,1))
data(wrld_simpl)
plot(wrld_simpl, add=TRUE,col='dark grey')

## PREDICTIVA RAC
plot(predm_glm_rac, col=tim.colors(100)[1:100],main="Glm RAC", axes=T)
data(wrld_simpl)
plot(wrld_simpl, add=TRUE,col='dark grey')


