library(terra)
library(tidyverse)

#Rare species
rare <- readRDS("rds/pa_rare.rds") %>% vect()
rare <- as.data.frame(rare)
rare <- rare[, 2:393] #392 raras
colnames(rare) <- gsub("\\.", " ", colnames(rare))
rare <- colnames(rare)

#Already filtered data from GBIF and CRIA
herb <- read_csv("data/tables/plantR_herbarium_filtered.csv")
herb <- herb[, c("scientificName.new1", 
                 "decimalLatitude.new1", 
                 "decimalLongitude.new1")]

rares_herb <- herb %>% filter(scientificName.new1 %in% rare)
n_distinct(rares_herb$scientificName.new1) #apenas 241 nos herbários
#Diferença possivelmente pela lista de coletores confiáveis
names(rares_herb) <- c("spp", "lat", "lon")
rares_herb <- rares_herb[, c(1, 3, 2)]
write.csv(rares_herb, "data/tables/herbarium_rares.csv")











