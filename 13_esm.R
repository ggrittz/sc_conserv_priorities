##### Modeling rare species #####
library(ecospat)
library(biomod2)
library(terra)
library(tidyverse)
library(predicts)

#Sites and environmental data
sites <- readRDS('rds/pa_methodENphylo.rds')
rares <- readRDS('data/pop_and_geo_rares.rds')$spp
sites <- as.data.frame(sites, geom = "XY")
sites_spp <- sites[, 1:81]
sites_fct <- sites[, 82:ncol(sites)]
sites_spp <- sites_spp %>% select(all_of(rares))
sites <- cbind(sites_spp, sites_fct)

env_terra <- readRDS('rds/predictors_raster.rds')

#List of species 
sp <- colnames(sites[, 1:41])

#Result data frame
#Data.frame to fill the results
result = structure(list(spp = character(), 
                        boyce = numeric(), 
                        class = "data.frame"))

#Species test using biomod2 and ecospat
for(i in sp) {
     
print(paste('Modeling', i), sep = ' ')

#Random points automatically avoid sampling bg on site points
bg <- backgroundSample(mask = env_terra[[1]],
                       n = 1000, #
                       p = sites[, c("x", "y")])

#Adding background points to site data
bg <- as.data.frame(bg)
#sites <- as.data.frame(sites, geom = "XY")
all_points <- bind_rows(bg, sites) #bg points plus sites

# Get corresponding presence/absence data
myResp <- as.numeric(all_points[, i])
# Transform true absences into potential pseudo-absences
myResp.PA <- ifelse(myResp == 1, 1, NA)

# Get corresponding XY coordinates
myRespXY <- all_points[, c('x', 'y')]

# Environmental variables (selecting only 3 out of 12 to model faster)
myExpl <- env_terra

# Format Data with pseudo-absences : random method
myBiomodData.r <- BIOMOD_FormatingData(resp.var = myResp.PA,
                                       expl.var = myExpl,
                                       resp.xy = myRespXY,
                                       resp.name = i,
                                       PA.nb.rep = 1,
                                       PA.nb.absences = 1000 + sum(is.na(myResp)), #selecting only 1000 now
                                       PA.strategy = 'random',
                                       filter.raster = TRUE)

#Default
bmOptions <- bm_DefaultModelingOptions()

# Modeling
my.ESM <- ecospat.ESM.Modeling(data = myBiomodData.r,
                               NbRunEval = 10,
                               DataSplit = 80,
                               weighting.score = "SomersD",
                               models = "ANN",
                               #tune = TRUE,
                               models.options = bmOptions,
                               parallel = TRUE)
                         
# Evaluation and average of simple bivariate models
my.ESM_EF <- ecospat.ESM.EnsembleModeling(ESM.modeling.output = my.ESM,
                                          weighting.score = "Boyce",
                                          threshold = 0.25)

evalBoyce <- mean(my.ESM_EF$ESM.evaluations$Boyce)

### Projection of simple bivariate models into new space
my.ESM_proj_current <- ecospat.ESM.Projection(ESM.modeling.output = my.ESM,
                                              new.env = myExpl)

### Projection of calibrated ESMs into new space
my.ESM_EFproj_current <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = my.ESM_proj_current,
                                                        ESM.EnsembleModeling.output = my.ESM_EF)

final_raster <- my.ESM_EFproj_current/1000
writeRaster(final_raster, 
            paste("results/esm/rasters/", i, ".tif", sep = ""),
            overwrite = TRUE)

#Evaluation table for each species
result = rbind(result,
               data.frame(spp = i,
                          boyce = evalBoyce))

mywd <- list.dirs(full.names = TRUE)
toremove <- mywd[grepl("ESM.BIOMOD.output", mywd)][1]
unlink(toremove, recursive = TRUE)

     }

write.csv(result, "results/esm/evaluation.csv")





