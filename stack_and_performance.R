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
#Richness maps
ensemble_model <- round(mean(rf_model, brt_model, gam_model, 
                             glm_model, mars_model, bart_model, 
                             svm_model, xgb_model), 0)

#SSDM means matrix
ensemble_model_matrix <- round(mean(rf_model_matrix, brt_model_matrix, gam_model_matrix,
                                    glm_model_matrix, mars_model_matrix, bart_model_matrix,
                                    svm_model_matrix, xgb_model_matrix), 2)

#Ensemble standard deviation
#Richness maps
ensemble_model_sd <- round(stdev(rf_model, brt_model, gam_model,
                                 glm_model, mars_model, bart_model, 
                                 svm_model, xgb_model), 0)

#SSDM means matrix
ensemble_model_matrix_sd <- round(stdev(rf_model_matrix, brt_model_matrix, gam_model_matrix,
                                        glm_model_matrix, mars_model_matrix, bart_model_matrix,
                                        svm_model_matrix, xgb_model_matrix), 2)


#Save matrices
writeRaster(rf_model_matrix, filename = "maps/rf_map_matrix.tif", overwrite = TRUE)
writeRaster(brt_model_matrix, filename = "maps/brt_map_matrix.tif", overwrite = TRUE)
writeRaster(gam_model_matrix, filename = "maps/gam_map_matrix.tif", overwrite = TRUE)
writeRaster(glm_model_matrix, filename = "maps/glm_map_matrix.tif", overwrite = TRUE)
writeRaster(mars_model_matrix, filename = "maps/mars_map_matrix.tif", overwrite = TRUE)
writeRaster(bart_model_matrix, filename = "maps/bart_map_matrix.tif", overwrite = TRUE)
writeRaster(svm_model_matrix, filename = "maps/svm_map_matrix.tif", overwrite = TRUE)
writeRaster(xgb_model_matrix, filename = "maps/xgb_map_matrix.tif", overwrite = TRUE)
writeRaster(ensemble_model_matrix_sd, filename = "maps/ensemble_model_matrix_sd.tif",
            overwrite = TRUE)
writeRaster(ensemble_model_matrix, filename = "maps/ensemble_model_matrix.tif",
            overwrite = TRUE)


#rm(list=ls())

#################################################
#################################################
#################################################
##### Ranking models by AUC-ROC performance #####
#################################################
#################################################
library(ggplot2)
library(tidyverse)

models <- c("bart", "brt", "gam", "glm", "mars", "rf", "svm", "xgboost")

r <- vector('list', length(models))

for (i in seq_along(models)){
     r[[i]] <- read.csv(paste0("results/", models[i], "/", "evaluation.csv"))
}

##### Adjusting names #####
bart_eval <- r[[1]]
names(bart_eval) <- c("model", "spp", "roc", "prc")
bart_eval$model <- "bart"

brt_eval <- r[[2]]
names(brt_eval) <- c("model", "spp", "roc", "prc")
brt_eval$model <- "brt"

gam_eval <- r[[3]]
names(gam_eval) <- c("model", "spp", "roc", "prc")
gam_eval$model <- "gam"

glm_eval <- r[[4]]
names(glm_eval) <- c("model", "spp", "roc", "prc")
glm_eval$model <- "glm"

mars_eval <- r[[5]]
names(mars_eval) <- c("model", "spp", "roc", "prc")
mars_eval$model <- "mars"

rf_eval <- r[[6]]
names(rf_eval) <- c("model", "spp", "roc", "prc")
rf_eval$model <- "rf"

svm_eval <- r[[7]]
names(svm_eval) <- c("model", "spp", "roc", "prc")
svm_eval$model <- "svm"

xgb_eval <- r[[8]]
names(svm_eval) <- c("model", "spp", "roc", "prc")
xgb_eval$model <- "xgb"
xgb_eval <- xgb_eval[, c(2:5)]


all_eval <- rbind.data.frame(bart_eval, brt_eval, gam_eval, glm_eval, 
                             mars_eval, rf_eval, svm_eval, xgb_eval)

all_eval <- all_eval[, c(1, 3, 4)]

#Grouping to extract the AUC-ROC and AUC-PRC mean and sd
group_summ <- all_eval %>%
     group_by(model) %>%
     summarise(roc_mean = round(mean(roc), 2),
               roc_sd = round(sd(roc), 2),
               prc_mean = round(mean(prc), 2),
               prc_sd = round(sd(prc), 2))

#Plotting (add median roc e median prc depois)
ggplot(all_eval, aes(x = roc, y = prc)) +
     geom_point(alpha = 0.5, size = 1) +
     #geom_point(data = group_summ, size = 3) +
     labs(x = "ROC", y = "PRC") +
     geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
     geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
     theme(legend.position = "bottom") +
     theme_bw()
     

ggsave(filename = "figureS1.svg",
       plot = last_plot(),
       device = "svg",
       path = "paper/Figuras/",
       width = 5,
       height = 3,
       bg = NULL)


##### PREDICTED x OBSERVED #####
#Observed data
sites <- readRDS('rds/pa_common.rds')
sites <- terra::as.data.frame(sites, xy = TRUE, geom = 'XY')
sites_new <- sites$ua
observed <- rowSums(sites[, 2:279])
sites <- vect(sites, geom = c("x", "y"))

#Predicted data
predicted_rf <- terra::extract(x = rf_model, y = sites)$sum
predicted_brt <- terra::extract(x = brt_model, y = sites)$sum
predicted_gam <- terra::extract(x = gam_model, y = sites)$sum
predicted_glm <- terra::extract(x = glm_model, y = sites)$sum
predicted_mars <- terra::extract(x = mars_model, y = sites)$sum
predicted_bart <- terra::extract(x = bart_model, y = sites)$sum
predicted_svm <- terra::extract(x = svm_model, y = sites)$sum
predicted_xgb <- terra::extract(x = xgb_model, y = sites)$sum
predicted_ensemble <- terra::extract(x = ensemble_model, y = sites)$sum

eval_df <- cbind.data.frame(predicted_rf, predicted_brt, predicted_gam,
                            predicted_glm, predicted_mars, predicted_bart,
                            predicted_svm, predicted_xgb, predicted_ensemble, 
                            observed)

names(eval_df) <- c("RF", "BRT", "GAM", "GLM", "MARS",
                    "BART", "SVM", "XGBoost", "Ensemble",
                    "Observed")

#Each plot
list_pred <- vector(mode = "list",
                    length = 9) #9 models
names(list_pred) <- names(eval_df[, 1:9])


#List of all regression plots
for (i in names(eval_df[, 1:9])) {
     message(i)
     list_pred[[i]] <- local({
          
          p1 <- ggplot(data = eval_df,
                       mapping = aes(x = observed, y = .data[[i]])) + 
               geom_point(size = 0.5) +
               geom_smooth(method='lm', formula = y ~ x) +
               ggpmisc::stat_poly_eq() + 
               theme(axis.text = element_text(size = 10),
                     axis.title = element_text(size = 14),
                     axis.title.x = element_blank()) +
               ylab(i)
          print(p1)
     })
}

ggpred <- annotate_figure(ggarrange(plotlist = list_pred), 
                bottom = text_grob("Observed", size = 14))

ggpred

ggsave(filename = "figureS2.pdf",
       plot = last_plot(),
       device = "pdf",
       path = "paper/Figuras/",
       width = 8,
       height = 6,
       bg = NULL)



