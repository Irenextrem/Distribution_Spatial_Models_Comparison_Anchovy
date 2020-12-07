# PROBLEMA: 

### PRESENCIA DE ANCHOAS ###
data <- read.csv('Anchoas_Aqua.csv',TRUE,",") 
colnames(data)<- c('Género','Especie','Lat','Lon','C.Square','Profundidad (metros)','Tª sup agua','Tª prof agua','Salinidad sup (psu)','Salinidad prof (psu)','Producción Primaria (mgC.m.3.day1)','Concentración de hielo en mar (0.1 fracción)','Oxígeno disuelto (mmol.m.3)', 'Distancia a la costa (km)')
data <- data[c(-151, -152, -153, -160, -163, -164, -165, -166),]
dups2 <- duplicated(data[, c("Especie","Lon","Lat")]) #Compruebo que no hay duplicados en esas columnas
sum(dups2) #Veo cuántos duplicados son: 258
data <- data[!dups2, ] #Verificar si son diferentes antes de eliminarlos.

### VARIABLES AMBIENTALES ###
files <-(list.files("D:/Desktop/Remember/Estudios/Educación Formal/Máster/Máster Valencia/Bioestadística/Curso 2/Especialización/8 Distribución de Especies/Clase 1/Practice1/predictors", full.names=T, pattern=".tif"))
predictors <- stack(files)
names(predictors) <- c("calcite","chlomean","nitrate","ph","phos","salinity","silicate","sstmean")

# Empiezo a tratar los datos
ext<-extent(-25,40,15,80)
predictors<-crop(predictors,ext)
plot(predictors)

# Estandarizo las variables
predictors2 <- scale(predictors)

### COORDENADAS ###
coords<- cbind(data$Lon,data$Lat)
colnames(coords)<- c('x','y')

### SDMDATA ###
backgr <- randomPoints(predictors2, 855)

presvals <- extract(predictors2, coords)  #Hay muchos NAs ¿Qué pasa?
absvals <- extract(predictors2, backgr)
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals))) #1 donde hay presencia y 0 donde hay ausencias

sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals), coords)) # NO ME DEJA CREARLO 

### RANDOM FOREST ###

model <- pb ~ calcite+nitrate+salinity
rf3 <- randomForest(model, sdmdata);rf3
varImpPlot(rf3)
