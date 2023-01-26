##############################################
##############################################
##### Stacking species distribution maps #####
##############################################
##############################################
library(terra)

#source("00_ecospat_prr.R")

models <- c("rf", "brt", "gam", "glm", "mars", "bart", "svm", "xgboost")

r <- vector('list', length(models))

for (i in seq_along(models)){
     r[[i]] <- list.files(paste0("results/", models[i], "/rasters/"), 
                          pattern="tif$", full.names=TRUE)
     }


#Sum of rough probabilities from s-SDM
rf_model <- round(sum(rast(r[[1]])), 0)
rf_model_matrix <- rast(r[[1]])

brt_model <- round(sum(rast(r[[2]])), 0)
brt_model_matrix <- rast(r[[2]])

gam_model <- round(sum(rast(r[[3]])), 0)
gam_model_matrix <- rast(r[[3]])

glm_model <- round(sum(rast(r[[4]])), 0)
glm_model_matrix <- rast(r[[4]])

mars_model <- round(sum(rast(r[[5]])), 0)
mars_model_matrix <- rast(r[[5]])

bart_model <- round(sum(rast(r[[6]])), 0)
bart_model_matrix <- rast(r[[6]])

svm_model <- round(sum(rast(r[[7]])), 0)
svm_model_matrix <- rast(r[[7]])

xgb_model <- round(sum(rast(r[[8]])), 0)
xgb_model_matrix <- rast(r[[8]])

#Ensemble (unweighted mean of all models)
ensemble_model <- round(mean(rf_model, brt_model, gam_model, 
                             glm_model, mars_model, bart_model, 
                             svm_model, xgb_model), 0)
#Ensemble standard deviation
ensemble_model_sd <- round(stdev(rf_model, brt_model, gam_model,
                                 glm_model, mars_model, bart_model, 
                                 svm_model, xgb_model), 0)

#Visualizing
#Random Forest
plot(rf_model, main = "Random Forest", cex = 2)
writeRaster(rf_model, filename = "maps/rf_map.tif", overwrite = TRUE)

#Boosting Regression Trees
plot(brt_model, main = "Boosting Regression Trees", cex = 2)
writeRaster(brt_model, filename = "maps/brt_map.tif", overwrite = TRUE)

#Generalized Additive Model
plot(gam_model, main = "Generalized Additive Model", cex = 2)
writeRaster(gam_model, filename = "maps/gam_map.tif", overwrite = TRUE)

#Generalized Linear Model
plot(glm_model, main = "Generalized Linear Model", cex = 2)
writeRaster(glm_model, filename = "maps/glm_map.tif", overwrite = TRUE)

#Multivariate Adaptative Regression Spline
plot(mars_model, main = "Multivariate Adaptative Regression Spline", cex = 2)
writeRaster(mars_model, filename = "maps/mars_map.tif", overwrite = TRUE)

#Bayesian Additive Regression Trees
plot(bart_model, main = "Bayesian Additive Regression Trees", cex = 2)
writeRaster(bart_model, filename = "maps/bart_map.tif", overwrite = TRUE)

#Support Vector Machine
plot(svm_model, main = "Support Vector Machine", cex = 2)
writeRaster(svm_model, filename = "maps/svm_map.tif", overwrite = TRUE)

#eXtreme Gradient Boosting
plot(xgb_model, main = "eXtreme Gradient Boosting", cex = 2)
writeRaster(xgb_model, filename = "maps/xgb_map.tif", overwrite = TRUE)

#Ensemble (unweighted mean)
plot(ensemble_model, main = "Ensemble Average", cex = 2)
writeRaster(ensemble_model, filename = "maps/ensemble_map.tif", overwrite = TRUE)

#Ensemble standard deviation
plot(ensemble_model_sd, main = "Ensemble Standard Deviation", cex = 2)
writeRaster(ensemble_model_sd, filename = "maps/ensemble_sd_map.tif", overwrite = TRUE)

#Save matrices
writeRaster(rf_model_matrix, filename = "maps/rf_map_matrix.tif", overwrite = TRUE)
writeRaster(brt_model_matrix, filename = "maps/brt_map_matrix.tif", overwrite = TRUE)
writeRaster(gam_model_matrix, filename = "maps/gam_map_matrix.tif", overwrite = TRUE)
writeRaster(glm_model_matrix, filename = "maps/glm_map_matrix.tif", overwrite = TRUE)
writeRaster(mars_model_matrix, filename = "maps/mars_map_matrix.tif", overwrite = TRUE)
writeRaster(bart_model_matrix, filename = "maps/bart_map_matrix.tif", overwrite = TRUE)
writeRaster(svm_model_matrix, filename = "maps/svm_map_matrix.tif", overwrite = TRUE)
writeRaster(xgb_model_matrix, filename = "maps/xgb_map_matrix.tif", overwrite = TRUE)


rm(list=ls())









