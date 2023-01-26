#######################################################
#######################################################
##### Modeling Bayesian Additive Regression Trees #####
#######################################################
#######################################################

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
          
          #Separate training-testing by pre-defined spatial blocks
          split <- initial_split(sites_sp, prop = .7, strata = block) #70/30 in each block
          train <- training(split)
          test <- testing(split)
          
          #Running the BART model
          #Model training here
          bart_train <- bart.step(x.data = train[, 5:9],
                                  y.data = train$occ)
          
          #Model testing here
          test_pred <- dbarts:::predict.bart(bart_train, test)
          
          #Since 1000 models are runned we need to extract the mean
          test_pred <- colMeans(test_pred)
          
          #ROC-AUC and PRC-AUC evaluation
          precrec_obj <- evalmod(scores = test_pred,
                                 labels = test$occ)
          
          #Prediction maps based on trained data
          bartpred <- predict2.bart(object = bart_train,
                              x.layers = env_raster,
                              splitby = 20,
                              quiet = TRUE)
          
          #Saving the output prediction raster
          writeRaster(bartpred,
                      paste("results/bart/rasters/", s, ".tif", sep = ""),
                      overwrite = TRUE)
          
          #Evaluation table for each species
          result = rbind(result,
                         data.frame(spp = s,
                                    roc = round(precrec::auc(precrec_obj)$aucs[1], 4),
                                    prc = round(precrec::auc(precrec_obj)$aucs[2], 4)))
          
          
     } else {
          
          
          #Separate training-testing by pre-defined spatial blocks
          split <- initial_split(sites_sp, prop = .7) #70/30 in each block
          train <- training(split)
          test <- testing(split)
          
          #Running the BART model
          #Model training here
          bart_train <- bart.step(x.data = train[, 5:9],
                                  y.data = train$occ)
          
          #Model testing here
          test_pred <- dbarts:::predict.bart(bart_train, test)
          
          #Since 1000 models are runned we need to extract the mean
          test_pred <- colMeans(test_pred)
          
          #ROC-AUC and PRC-AUC evaluation
          precrec_obj <- evalmod(scores = test_pred,
                                 labels = test$occ)
          
          #Prediction maps based on trained data
          bartpred <- predict2.bart(object = bart_train,
                                    x.layers = env_raster,
                                    splitby = 20,
                                    quiet = TRUE)
          
          #Saving the output prediction raster
          writeRaster(bartpred,
                      paste("results/bart/rasters/", s, ".tif", sep = ""),
                      overwrite = TRUE)
          
          #Evaluation table for each species
          result = rbind(result,
                         data.frame(spp = s,
                                    roc = round(precrec::auc(precrec_obj)$aucs[1], 4),
                                    prc = round(precrec::auc(precrec_obj)$aucs[2], 4)))
          
     }             
     
}

#Table
write.csv(result, "results/bart/evaluation.csv")

}