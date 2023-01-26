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
xgb_eval$model <- "xbg"
xgb_eval <- xgb_eval[, c(2:5)]


all_eval <- rbind.data.frame(bart_eval, brt_eval, gam_eval, glm_eval, 
                             mars_eval, rf_eval, svm_eval, xgb_eval)

all_eval <- all_eval[, c(1, 3, 4)]

#Grouping to extract the AUC-ROC and AUC-PRC mean
group_mean <- all_eval %>%
     group_by(model) %>%
     summarise(roc = mean(roc),
               prc = mean(prc))

#Plotting (add median roc e median prc depois)
ggplot(all_eval, aes(x = roc, y = prc, color = model)) +
     geom_point(alpha = 0.2) +
     geom_point(data = group_mean, size = 3) +
     theme_bw()


##### PREDICTED x OBSERVED #####
source("00_loading_data.R")

sites_new <- sites$ua
observed <- rowSums(sites[, 2:279])
sites <- vect(sites, geom = c("x", "y"))
predicted_rf <- terra::extract(x = rf_model, y = sites)$sum
predicted_brt <- terra::extract(x = brt_model, y = sites)$sum
predicted_gam <- terra::extract(x = gam_model, y = sites)$sum
predicted_glm <- terra::extract(x = glm_model, y = sites)$sum
predicted_mars <- terra::extract(x = mars_model, y = sites)$sum
predicted_bart <- terra::extract(x = bart_model, y = sites)$sum
predicted_svm <- terra::extract(x = svm_model, y = sites)$sum
predicted_xgb <- terra::extract(x = xgb_model, y = sites)$sum

eval_df <- cbind.data.frame(predicted_rf, predicted_brt, predicted_gam,
                            predicted_glm, predicted_mars, predicted_bart,
                            predicted_svm, observed)

svm_pred <- ggplot(data = eval_df,
                   mapping = aes(x = observed, y = predicted_svm)) + geom_point() +
          geom_smooth(method='lm', formula = y ~ x) +
     theme(axis.text = element_text(size = 12),
           axis.title = element_text(size = 16))




