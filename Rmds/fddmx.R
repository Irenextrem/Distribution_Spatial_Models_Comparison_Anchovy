# Función definitiva definitivísima milenaria Xtrem para hacer una cross validation 10 veces y que te devuelva
# al módico precio de un par de segundos una matriz con: AUC, COR, kappa, Sensitivity, specificity
# y el TSS

fddmx <- function(coordenadas, predictores,background){
  
  group <- kfold(coordenadas, 5)
  
  pres_train <- coordenadas[group != 1, ] #Las que no sean 1
  pres_test <- coordenadas[group == 1, ] #Las que sean 1
  
  bc <- bioclim(predictores, pres_train)
  
  group <- kfold(background, 5) #Hago 5 grupos de esos puntos
  backg_train <- background[group != 1, ]
  backg_test <- background[group == 1, ]
  
  eval.modesta <- evaluate(pres_test, backg_test, bc, predictores)
  
  auc_bc.model <- eval.modesta@auc #auc
  
  cor_bc.model <- eval.modesta@cor #cor
  
  kappa_bc.model <- mean(eval.modesta@kappa) #Kappa
  
  sensibility_bc.model <- mean(eval.modesta@TPR/(eval.modesta@TPR+eval.modesta@FNR)) #Sensibilidad
  
  specificity_bc.modelo <- mean(eval.modesta@TNR/(eval.modesta@FPR+eval.modesta@TNR)) #Especificidad
  
  TSSbc.model <- mean(eval.modesta@TPR+eval.modesta@TNR-1) #TSS
  
  return(c(auc_bc.model,cor_bc.model,kappa_bc.model,sensibility_bc.model,specificity_bc.modelo,TSSbc.model))
}