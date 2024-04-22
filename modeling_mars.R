##############################################################
##############################################################
##### Modeling Multivariate Adaptative Regression Spline #####
##############################################################
##############################################################

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
          
          #Response must be a factor
          train$occ <- as.factor(train$occ)
          levels(train$occ) <- c("C0", "C1")
          
          #Creating a tuning grid
          tune_grid <- expand.grid(degree = 1:3,
                                   nprune = seq(2, 100, length.out = 10) %>% floor()
                                   )
          
          #Using caret for cross-validation
          train_control <- trainControl(method = "cv",
                                        number = 10,
                                        classProbs = TRUE,
                                        summaryFunction = twoClassSummary,
                                        allowParallel = TRUE)
          
          #Starting cluster
          cluster <- makeCluster(8)
          registerDoParallel(cluster)
          
          #Running the MARS model
          #Model training here
          mars_fit <- caret::train(x = train[, 5:9],
                                   y = train$occ,
                                   method = "earth",
                                   metric = "ROC",
                                   trControl = train_control,
                                   tuneGrid = tune_grid)
          
          #Ending cluster
          stopCluster(cluster)
          registerDoSEQ()
          
          #Model testing here
          test_pred <- predict(mars_fit,
                               test,
                               type = "prob")$C1
          
          #ROC-AUC and PRC-AUC evaluation
          precrec_obj <- evalmod(scores = test_pred,
                                 labels = test$occ)
          
          #Prediction maps based on trained data
          marspred <- raster::predict(model = mars_fit,
                                     object = env_raster,
                                     type = "prob",
                                     index = 2)
          
          #Saving the output prediction raster
          writeRaster(marspred,
                      paste("results/mars/rasters/", s, ".tif", sep = ""),
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
          
          #Response must be a factor for {caret}
          train$occ <- as.factor(train$occ)
          levels(train$occ) <- c("C0", "C1")
          
          #Creating a tuning grid
          tune_grid <- expand.grid(degree = 1:3,
                                   nprune = seq(2, 100, length.out = 10) %>% 
                                        floor())
          
          #Using caret for cross-validation
          train_control <- trainControl(method = "cv",
                                        number = 10,
                                        classProbs = TRUE,
                                        summaryFunction = twoClassSummary,
                                        allowParallel = TRUE)
          
          #Starting cluster
          cluster <- makeCluster(8)
          registerDoParallel(cluster)
          
          #Running the MARS model
          #Model training here
          mars_fit <- caret::train(x = train[, 5:9],
                                   y = train$occ,
                                   method = "earth",
                                   metric = "ROC",
                                   trControl = train_control,
                                   tuneGrid = tune_grid)
          
          #Ending cluster
          stopCluster(cluster)
          registerDoSEQ()
          
          #Model testing here
          test_pred <- predict(mars_fit,
                               test,
                               type = "prob")$C1
          
          #ROC-AUC and PRC-AUC evaluation
          precrec_obj <- evalmod(scores = test_pred,
                                 labels = test$occ)
          
          #Prediction maps based on trained data
          marspred <- raster::predict(model = mars_fit,
                                      object = env_raster,
                                      type = "prob",
                                      index = 2)
          
          #Saving the output prediction raster
          writeRaster(marspred,
                      paste("results/mars/rasters/", s, ".tif", sep = ""),
                      overwrite = TRUE)
          
          #Evaluation table
          result = rbind(result,
                         data.frame(spp = s,
                                    roc = round(precrec::auc(precrec_obj)$aucs[1], 4),
                                    prc = round(precrec::auc(precrec_obj)$aucs[2], 4)))
          
     }
     
}

#Table
write.csv(result, "results/mars/evaluation.csv")

}


