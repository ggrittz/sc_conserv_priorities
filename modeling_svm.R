###########################################
###########################################
##### Modeling Support Vector Machine #####
###########################################
###########################################

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
          
          #Defining class weights
          pr_num <- as.numeric(table(train$occ)["C1"]) #N of presences
          ab_num <- as.numeric(table(train$occ)["C0"]) #N of absences
          cwt <- c("C0" = pr_num / ab_num, "C1" = 1)
          
          #Creating a tuning grid
          tune_grid <- expand.grid(C = c(0.01, 0.1, 1, 10, 100),
                                   sigma = c(0.001, 0.01, 0.1))
          
          #Using caret for cross-validation
          train_control <- trainControl(method = "cv",
                                        number = 10,
                                        classProbs = TRUE,
                                        summaryFunction = twoClassSummary,
                                        allowParallel = TRUE)
          
          #Starting cluster
          cluster <- makeCluster(8)
          registerDoParallel(cluster)
          
          #Running the SVM model
          #Model training here
          svm_fit <- caret::train(x = train[, 5:9],
                                  y = train$occ,
                                  method = "svmRadialSigma",
                                  metric = "ROC",
                                  trControl = train_control,
                                  tuneGrid = tune_grid,
                                  class.weights = cwt)
          
          #Ending cluster
          stopCluster(cluster)
          registerDoSEQ()
          
          #Model testing here
          test_pred <- predict(object = svm_fit,
                               newdata = test[, 5:9],
                               type = "prob")$C1
          
          #ROC-AUC and PRC-AUC evaluation
          precrec_obj <- evalmod(scores = test_pred,
                                 labels = test$occ)
          
          #Prediction maps based on trained data
          svmpred <- raster::predict(model = svm_fit,
                                      object = env_raster,
                                      type = "prob",
                                     index = 2)
          
          #Saving the output prediction raster
          writeRaster(svmpred,
                      paste("results/svm/rasters/", s, ".tif", sep = ""),
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
          
          #Response must be a factor
          train$occ <- as.factor(train$occ)
          levels(train$occ) <- c("C0", "C1")
          
          #Defining class weights
          pr_num <- as.numeric(table(train$occ)['C1']) #N of presences
          ab_num <- as.numeric(table(train$occ)['C0']) #N of absences
          cwt <- c("C0" = pr_num / ab_num, "C1" = 1)
          
          #Creating a tuning grid
          tune_grid <- expand.grid(C = c(0.01, 0.1, 1, 10, 100),
                                   sigma = c(0.001, 0.01, 0.1))
          
          #Using caret for cross-validation
          train_control <- trainControl(method = "cv",
                                        number = 10,
                                        classProbs = TRUE,
                                        summaryFunction = twoClassSummary,
                                        allowParallel = TRUE)
          
          #Starting cluster
          cluster <- makeCluster(8)
          registerDoParallel(cluster)
          
          #Running the SVM model
          #Model training here
          svm_fit <- caret::train(x = train[, 5:9],
                                  y = train$occ,
                                  method = "svmRadialSigma",
                                  metric = "ROC",
                                  trControl = train_control,
                                  tuneGrid = tune_grid,
                                  class.weights = cwt)
          
          #Ending cluster
          stopCluster(cluster)
          registerDoSEQ()
          
          #Model testing here
          test_pred <- predict(object = svm_fit,
                               newdata = test[, 5:9],
                               type = "prob")$C1
          
          #ROC-AUC and PRC-AUC evaluation
          precrec_obj <- evalmod(scores = test_pred,
                                 labels = test$occ)
          
          #Prediction maps based on trained data
          svmpred <- raster::predict(model = svm_fit,
                                     object = env_raster,
                                     type = "prob",
                                     index = 2)
          
          #Saving the output prediction raster
          writeRaster(svmpred,
                      paste("results/svm/rasters/", s, ".tif", sep = ""),
                      overwrite = TRUE)
          
          #Evaluation table for each species
          result = rbind(result,
                         data.frame(spp = s,
                                    roc = round(precrec::auc(precrec_obj)$aucs[1], 4),
                                    prc = round(precrec::auc(precrec_obj)$aucs[2], 4)))
          
     }
     
}

#Evaluation
write.csv(result, 'results/svm/evaluation.csv')

}
