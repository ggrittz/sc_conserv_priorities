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
sites <- readRDS('rds/pa_common.rds')
sites <- terra::unwrap(sites)

#This is a site-by-species matrix of 0s and 1s
sites <- terra::as.data.frame(sites, xy = TRUE, geom = 'XY')

#Raster cells over the study area (Santa Catarina)
env_terra <- readRDS('rds/predictors_df.rds') %>% rast()
env_raster <- brick(env_terra) #for blockCV
#env <- terra::as.data.frame(env, xy = TRUE, geom = 'XY')

#List of species
sp <- colnames(sites[, 2:279])

#Data.frame to fill the results
result = structure(list(spp = character(), 
                        roc = numeric(), 
                        prc = numeric()), 
                   class = "data.frame")
