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

# Base de datos sin escalar
sdmdata_nos <- read.csv("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Datos/sdmdata_sin_escalar.csv") #Sin las variables escaladas
sdmdata_nos<- sdmdata_nos[,c(-1)]

#GAM sin escalar
gam1 <- gam(pb ~ s(salinity)+s(bathy)+s(tempmean)+s(odismean), family = binomial, data =sdmdata_nos) 

#Predictiva

pred_gam<-predict(predictors3,gam1,type='response')

saveRDS(pred_gam, "C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/GAM/Predicciones_gam_Atl_sinescalar.ascii")

#Representación

ggplot() +
  geom_raster(data = raster::as.data.frame(pred_gam , xy = TRUE) , aes(x = x, y = y, fill = layer)) +
  coord_equal() +
  labs( x = "", y = "")+theme_minimal()+
  scale_color_brewer(palette = 'YlGnBu')+ labs(fill='')

plot(pred_gam, col=tim.colors(100)[1:100],main="GAM", axes=T)
data(wrld_simpl)
plot(wrld_simpl, add=TRUE,col='dark grey')

########    RAC   ###################
#Extract residuals from the GAM model and map them
# r<- raster(xmn=-11, xmx=9, ymn=43, ymx=60, nrows=324, ncols=660);r #Hago un raster grande
# res(r) <‐ 0.083 #Le doy una resolución de 0.083
# xy <-cbind(sdmdata_nos$x, sdmdata_nos$y)#Cojo las coordenadas
# xy_residuals <-cbind(xy, resid(gam1)) #Uno las coordenadas y los residuos a un mismo objeto
# # predictors3
# 
# par(mfrow=c(1,2))
# r[cellFromXY(r, xy_residuals)] <-xy_residuals[,3] #Doy esos valoes de residuos a esas coordenadas
# plot(r,col='red') #Efectivamente me salen residuos
# ext <- c(-11,9,43,60) #Extensión a cortar
# salinity <- raster("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Predictores/RAC/GAM/salinity_atln.TIF") #Cargo este raster para poder ajustar el raster r
# r <- resample(r,salinity) #Ajusto las dimensiones de r a las de salinity
# r<-crop(r,ext);r #Lo corto
# # writeRaster(r, filename="RAC_gam_atln.tif", format="GTiff", overwrite=TRUE) #Lo guado en el wd
# 
# #Calculate residuals autocovariate
# focal_rac_rast <-focal(r, w=matrix(1,3,3), fun = mean,  na.rm = TRUE)
# 
# #Extract the values of the focal operation from focal_rac_rest raserfile using the coordinates stored in xy
# focal_rac_vect <-extract(focal_rac_rast, xy)
# length(focal_rac_vect)
# 
# #Add as a column to the data
# dd<-cbind(sdmdata_nos, focal_rac_vect)
# I <- is.na(dd$focal_rac_vect)
# ddd<- dd[!I,]
# dim(dd)
# dim(ddd)
# 
# write.csv(ddd, 'C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Datos/rac_gam_sinescalar.csv')
ddd<- read.csv('C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Datos/rac_gam_sinescalar.csv')
##############################
### --- Kriging con RAC--- ###
##############################

# salinity<- load_layers(c("BO_salinity"))# datadir = tempdir())
# e=extent(-11,9,43,60) #Extensión con la que trabajo
# sal=crop(salinity,e) #Corto el raster
# 
# # Lo pinto
# plot(sal)
# points(ddd$x,ddd$y)
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
# matrix<-cbind(ddd$focal_rac_vect,coords) #Uno las coordenadas a los valores de rac
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
# rast<- raster(xmn=-11, xmx=9, ymn=43, ymx=60, nrows=276, ncols=624)#2  
# p<- rasterize(xy, rast, new$depth, fun=max,na.rm=F) #Hago un raster
# p<-resample(p,sal) #Le doy las dimensiones de salinidad
# e<-extent(p) #
# 
# sp<-crop(RAC,e)
# sp=resample(sp,p)
# sp<-raster::mask(sp,p)


# Plot
# plot(sp,col=tim.colors(100)[1:100],main=" ", axes=T)
# data(wrld_simpl)
# plot(wrld_simpl, add=TRUE,col='dark grey')
# 
# plot(RAC,col=tim.colors(100)[1:100],main=" ", axes=T)
# data(wrld_simpl)
# plot(wrld_simpl, add=TRUE,col='dark grey')
# 
# writeRaster(sp, filename="sp_gam_atln_sinesc.tif", format="GTiff", overwrite=TRUE)

#Predictores
files<-(list.files("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Predictores/RAC/GAM", full.names=T, pattern=".tif"))
predictors <- stack(files)
names(predictors) <- c("bathy","odismean","salinity","focal_rac_vect","tempmean")
predictors_rac <- scale(predictors) 

#GAM RAC sin escalar
gam_rac <- gam(pb ~ s(salinity)+s(bathy)+s(tempmean)+s(odismean) + s(focal_rac_vect), family = binomial, data =ddd)

#Predictiva

pred_gam_rac <-predict(predictors_rac, gam_rac,type='response')
saveRDS(pred_gam_rac, "C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/GAM/Predicciones_gam_rac_Atl_sinescalar.ascii")

#Representación

ggplot() +
  geom_raster(data = raster::as.data.frame(pred_gam_rac , xy = TRUE) , aes(x = x, y = y, fill = layer)) +
  coord_equal() +
  labs( x = "", y = "")+theme_minimal()+
  scale_color_brewer(palette = 'YlGnBu')+ labs(fill='') 

par(mfrow=c(1,1))
plot(pred_gam_rac, col=tim.colors(100)[1:100],main="GAM RAC", axes=T)
data(wrld_simpl)
plot(wrld_simpl, add=TRUE,col='dark grey')

#Gam bivariante sin escalar
gam_esp_bis_w <- gam(pb ~ s(salinity)+s(bathy)+s(tempmean)+s(odismean) + s(y,x,bs="gp",k=100,m=c(1,1)), family = binomial, data =sdmdata_nos) #Asi lo pone en WOOD

# Creo los vectores sobre los que voy a poner la predicción
# xx <-seq(-11,9,by = 0.1)
# yy <- seq(43,60,by=0.1)
# length(xx) #201 puntos
# 201^2 #40401
# 
# #Realizo las repeticiones
# x_rep <- rep(xx, 40401)
# y_rep <- rep(yy, 40401)
# length(y_rep)#6908571
# length(x_rep)
# 
# #Ordeno los valores
# y_rep<- matrix(y_rep,ncol=1,nrow=6908571)
# x_rep<- matrix(x_rep,ncol=1,nrow=6908571)
# 
# y_rep_or <- y_rep[order(y_rep)]
# 
# summary(prueba)
# #Genero una base de datos
# new_coord<-cbind(x_rep,y_rep_or)
# colnames(new_coord)<-c("x","y")
# 
# # Extraigo los valores de los predictores para esos valores
# bd_new <- extract(predictors3, new_coord)
# 
# #Lo junto todo en la nueva base de datos
# sdmdata_new <- data.frame(cbind(new_coord,bd_new))
# summary(sdmdata_new)
# 
# # Quito los NAs
# to.remove <- which(!complete.cases(sdmdata_new))
# sdmdata_new <- sdmdata_new[-to.remove,]
# summary(sdmdata_new)
# 
# write.csv(sdmdata_new,"C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Datos/sdmdata_gam_biv_pred_sinescalar.csv")
sdmdata_new<- read.csv("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Datos/sdmdata_gam_biv_pred_sinescalar.csv")

#Predicción
# pred_gam_bi_w_2 <- predict(gam_esp_bis_w, sdmdata_new, se.fit = T, type='response') 
# saveRDS(pred_gam_bi_w_2, "C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/GAM/Predicciones_gam_bi_new_data_sinescalar.ascii")
pred_gam_bi_w_2<- readRDS("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/GAM/Predicciones_gam_bi_new_data_sinescalar.ascii")

#Hago un raster
r<- raster(xmn=-11, xmx=9, ymn=43, ymx=60, nrows=324, ncols=660);r #Hago un raster grande
res(r) <‐ 0.083 #Le doy una resolución de 0.083
xy <-cbind(sdmdata_new$x, sdmdata_new$y)#Cojo las coordenadas
xy_residuals <-cbind(xy,pred_gam_bi_w_2$fit)  #Uno las coordenadas y los residuos a un mismo objeto
# predictors3

par(mfrow=c(1,2))
r[cellFromXY(r, xy_residuals)] <-xy_residuals[,3] #Doy esos valoes de residuos a esas coordenadas
plot(r,col='red') #Efectivamente me salen residuos
ext <- c(-11,9,43,60) #Extensión a cortar
salinity <- raster("C:/Users/Irene/Source/Repositorios/TFM-Irene-Extremera/Predictores/RAC/GAM/salinity_atln.TIF") #Cargo este raster para poder ajustar el raster r
r <- resample(r,salinity) #Ajusto las dimensiones de r a las de salinity
r<-crop(r,ext);r #Lo corto
rr<- aggregate(r,2)
writeRaster(rr, filename="gam_biv_sinescalar_atln.tif", format="GTiff", overwrite=TRUE)

# Mapas predictivos gam

par(mfrow=c(1,2))

plot(pred_gam, col=tim.colors(100)[1:100],main="GAM", axes=T,zlim=c(0,1))
data(wrld_simpl)
plot(wrld_simpl, add=TRUE,col='dark grey')

plot(rr, col=tim.colors(100)[1:100],main="GAM bivariante", axes=T,zlim=c(0,1))
data(wrld_simpl)
plot(wrld_simpl, add=TRUE,col='dark grey')

plot(pred_gam_rac, col=tim.colors(100)[1:100],main="GAM RAC", axes=T)
data(wrld_simpl)
plot(wrld_simpl, add=TRUE,col='dark grey')

