##### Loading packages and data #####

#Use library(sessionInfo) for getting packages version to submit

{
     library(randomForest)
     library(terra)
     library(precrec)
     library(blockCV)
     library(tidyverse)
     library(raster)
     library(sf)
     library(rsample)
     library(dismo)
     library(mgcv)
     library(maxnet)
     library(readxl)
     library(caret)
     library(earth)
     library(doParallel)
     library(gam)
     library(MASS)
     library(embarcadero)
}

#Sites and environmental data
sites <- readRDS('rds/pa_common.rds') %>% vect()
sites <- terra::as.data.frame(sites, xy = TRUE, geom = 'XY')
env_terra <- readRDS('rds/predictors_raster.rds') %>% rast()
env_raster <- brick(env_terra) #for blockCV
#env <- terra::as.data.frame(env, xy = TRUE, geom = 'XY')

#List of species
sp <- colnames(sites[, 2:279])

#Data.frame to fill the results
result = structure(list(spp = character(), 
                        roc = numeric(), 
                        prc = numeric()), 
                   class = "data.frame")