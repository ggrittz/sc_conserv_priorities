################################################
################################################
##### Modeling Generalized Additive Models #####
################################################
################################################

#Loading everything needed for modeling
source("00_loading_data.R")

{
     #Only modeling common species
     for(s in sp){

     
     print(paste('Modeling', s))
          
     sites_sp <- sites[, c("x", "y", "ua", s,
                           "PC1", "PC2", "PC3", "PC4", "PC5")]
     names(sites_sp)[names(sites_sp) == s] <- "occ"
     

          if(sum(sites_sp$occ) >= 30){
          
          
          pa_data <- st_as_sf(sites_sp,
                              coords = c("x", "y"),
                              crs = crs(env_raster))
          
          #Separating the data into blocks of 165 km size (no more spatial autocorrelation)
          sb <- spatialBlock(speciesData = pa_data,
                             species = "occ",
                             rasterLayer = env_raster,
                             theRange = 165000, #Size of the blocks (m) from prior SAC analysis
                             k = 5, #Arbitrary
                             selection = "random",
                             iteration = 200, #Find evenly dispersed folds
                             biomod2Format = FALSE,
                             verbose = FALSE,
                             progress = FALSE)
          
          #Adding the partitioned blocks info to the initial data set
          sites_sp$block <- sb$foldID
          
          #Separate training-testing by pre-defined spatial blocks
          split <- initial_split(sites_sp, prop = .7, strata = block) #70/30 in each block
          train <- training(split)
          test <- testing(split)
          
          #Calculating case weights
          pr_num <- as.numeric(table(train$occ)["1"]) #N of presences
          ab_num <- as.numeric(table(train$occ)["0"]) #N of absences
          wt <- ifelse(train$occ == 1, 1, pr_num / ab_num)
     
          #Defining smooth parameters in GAM
          form <- occ ~ s(PC1) + s(PC2) + s(PC3) + s(PC4) + s(PC5)

          #Fitting the GAM model
          gm <- mgcv::gam(formula = as.formula(form),
                          data = train,
                          family = binomial(link = "logit"),
                          weights = wt,
                          method = "REML")

          #Model testing here
          test_pred <- predict(gm,
                               test,
                               type = "response")

          #ROC-AUC and PRC-AUC evaluation
          precrec_obj <- evalmod(scores = test_pred,
                                 labels = test$occ)

          #Prediction maps based on trained data
          gmpred <- raster::predict(model = gm,
                                    object = env_raster,
                                    type = "response")
          
          #Saving the output prediction raster
          writeRaster(gmpred,
                      paste("results/gam/rasters/", s, ".tif", sep = ""),
                      overwrite = TRUE)
          
          #Evaluation table for each species
          result = rbind(result,
                    data.frame(spp = s,
                               roc = round(precrec::auc(precrec_obj)$aucs[1], 4),
                               prc = round(precrec::auc(precrec_obj)$aucs[2], 4)))


     } else {
          

     #Separate training-testing by pre-defined spatial blocks
     split <- initial_split(sites_sp, prop = .70) #70/30 in each block
     train <- training(split)
     test <- testing(split)
          
     #Calculating case weights
     pr_num <- as.numeric(table(train$occ)["1"]) #N of presences
     ab_num <- as.numeric(table(train$occ)["0"]) #N of absences
     wt <- ifelse(train$occ == 1, 1, pr_num / ab_num)
          
     #Defining smooth parameters in GAM
     form <- occ ~ s(PC1) + s(PC2) + s(PC3) + s(PC4) + s(PC5)
     
     #Fitting the GAM model
     gm <- mgcv::gam(formula = as.formula(form),
                     data = train,
                     family = binomial(link = "logit"),
                     weights = wt,
                     method = "REML")
          
     #Model testing here
     test_pred <- predict(gm,
                          test, 
                          type = "response")
          
     #ROC-AUC and PRC-AUC evaluation
     precrec_obj <- evalmod(scores = test_pred,
                            labels = test$occ)
          
     #Prediction maps based on trained data
     gmpred <- raster::predict(model = gm,
                               object = env_raster,
                               type = "response")
          
     #Saving the output prediction raster
     writeRaster(gmpred,
                 paste("results/gam/rasters/", s, ".tif", sep = ""),
                 overwrite = TRUE)
          
     #Evaluation table for each species
     result = rbind(result,
                    data.frame(spp = s,
                               roc = round(precrec::auc(precrec_obj)$aucs[1], 4),
                               prc = round(precrec::auc(precrec_obj)$aucs[2], 4)))
     }
              
}


#Table
write.csv(result, "results/gam/evaluation.csv")

}