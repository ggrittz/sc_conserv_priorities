##################################
##################################
##### Modeling Random Forest #####
##################################
##################################

#Loading everything needed for modeling
source("00_loading_data.R")

{
     #Only modeling "common" species
     for(s in sp){
     
          
     print(paste('Modeling', s))
     
     sites_sp <- sites[, c("x", "y", "ua", s,
                           "PC1", "PC2", "PC3", "PC4", "PC5")]
     names(sites_sp)[names(sites_sp) == s] <- "occ"


          if(sum(sites_sp$occ) >= 30){
          
          
          pa_data <- st_as_sf(sites_sp,
                              coords = c("x", "y"),
                              crs = crs(env_terra))

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

          #Convert the response to factor for producing class relative likelihood
          #Random forest only classifies factors
          sites_sp$occ <- as.factor(sites_sp$occ)

          #Separate training-testing by pre-defined spatial blocks
          split <- initial_split(sites_sp, prop = .7, strata = block) #70/30 in each block
          train <- training(split)
          test <- testing(split)
          
          #Calculating class weights 
          pr_num <- as.numeric(table(train$occ)["1"]) #N of presences
          ab_num <- as.numeric(table(train$occ)["0"]) #N of absences
          smpsize <- c("1" = pr_num, "0" = ab_num)

          #Running the RF model
          #Model training here
          rf_downsample <- randomForest(formula = occ ~.,
                                        data = train[, 4:9],
                                        ntree = 1000,
                                        sampsize = smpsize,
                                        replace = TRUE)

          #Model testing here
          test_pred <- predict(rf_downsample,
                               test,
                               type = "prob",)[, "1"] #only presences are needed

          #ROC-AUC and PRC-AUC evaluation
          precrec_obj <- evalmod(scores = test_pred,
                                 labels = test$occ)

          #Prediction maps based on trained data
          rfpred <- raster::predict(model = rf_downsample,
                                    object = env_raster,
                                    type = "prob",
                                    index = 2)

          #Saving the output prediction raster
          writeRaster(rfpred,
                      paste("results/random_forest/rasters/", s, ".tif", sep = ""),
                      overwrite = TRUE)

          #Evaluation table for each species
          result = rbind(result,
                         data.frame(spp = s,
                                    roc = round(precrec::auc(precrec_obj)$aucs[1], 4),
                                    prc = round(precrec::auc(precrec_obj)$aucs[2], 4)))
     

          } else {
          
          
               #Convert the response to factor for producing class relative likelihood
               #Random forest only classifies factors
               sites_sp$occ <- as.factor(sites_sp$occ)
               
               #Separate training-testing by pre-defined spatial blocks
               split <- initial_split(sites_sp, prop = .7) #70/30 in each partition
               train <- training(split)
               test <- testing(split)
               
               #Calculating class weights 
               pr_num <- as.numeric(table(train$occ)["1"]) #N of presences
               ab_num <- as.numeric(table(train$occ)["0"]) #N of absences
               smpsize <- c("1" = pr_num, "0" = ab_num)
               
               #Running the RF model
               #Model training here
               rf_downsample <- randomForest(formula = occ ~.,
                                             data = train[, 4:9],
                                             ntree = 1000,
                                             sampsize = smpsize,
                                             replace = TRUE)
               
               #Model testing here
               test_pred <- predict(rf_downsample,
                                    test,
                                    type = "prob",)[, "1"] #only presences are needed
               
               #ROC-AUC and PRC-AUC evaluation
               precrec_obj <- evalmod(scores = test_pred,
                                      labels = test$occ)
               
               #Prediction maps based on trained data
               rfpred <- raster::predict(model = rf_downsample,
                                         object = env_raster,
                                         type = "prob",
                                         index = 2)
               
               #Saving the output prediction raster
               writeRaster(rfpred,
                           paste("results/random_forest/rasters/", s, ".tif", sep = ""),
                           overwrite = TRUE)
               
               #Evaluation table for each species
               result = rbind(result,
                              data.frame(spp = s,
                                         roc = round(precrec::auc(precrec_obj)$aucs[1], 4),
                                         prc = round(precrec::auc(precrec_obj)$aucs[2], 4)))

     }             

}


#Table
write.csv(result, "results/random_forest/evaluation.csv")

}