##### Modeling rare species (presences between 10 (?) and 19 points) #####

source("00_ecospat_spatstat.R")

{
library(terra)
library(tidyverse)
library(raster)
library(sf)
library(RRdtn)
library(V.PhyloMaker2)
library(flora)
library(RRphylo)
library(dismo)
library(ecospat)
library(spatstat)
library(predicts)
}


#Loading survey data for both rares and common species
spp_rares <- terra::as.data.frame(readRDS("rds/pa_rare_enfa.rds"), 
                                  geom = "XY")

spp_common <- terra::as.data.frame(readRDS("rds/pa_common_enfa.rds"),
                                   geom = "XY")

#Loading predictors' extent/mask
preds <- readRDS("rds/predictors_df.rds")
preds <- rast(preds)
crs(preds) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

spp_all <- merge(spp_common, spp_rares[, c(1:393)], by = "ua")
spp_all <- spp_all[, order(names(spp_all))]
spp_all <- spp_all %>% relocate(c("ua", "x", "y", "PC1", "PC2", "PC3", "PC4", "PC5"))
#saveRDS(spp_all, "rds/all_species_enfa.rds")

#Select only species with 10+ occurrences
only_spp <- spp_all[, c(9:678)]
only_spp <- only_spp[, colSums(only_spp) >= 10]
spp_all <- cbind(spp_all[, c(1:8)], only_spp)

#Getting families from Flora do Brasil
all_flora <- names(spp_all[, c(9:368)])
all_flora <- gsub("\\.", " ", all_flora)
all_flora <- get.taxa(all_flora)

#Data.frame for V.PhyloMaker2
all_df <- data.frame(species = all_flora$original.search,
                     genus = sapply(strsplit(all_flora$original.search, " "), `[`, 1),
                     family = all_flora$family,
                     species.relative = NA,
                     genus.relative = NA,
                     row.names = NULL)

#Remove species with less than 10 occurrences
prev_tree <- ape::read.tree("sample.tre")
to_remove <- prev_tree$tip.label
to_remove <- gsub("_", " ", to_remove)

#Remove NA families
all_df <- all_df %>% drop_na(family)

#Remove species that are not going into the tree (previously performed)
all_df <- all_df %>% filter(species %in% to_remove)

#V.PhyloMaker2
all_tree <- phylo.maker(sp.list = all_df,
                        tree = GBOTB.extended.TPL,
                        nodes = nodes.info.1.TPL,
                        scenarios = "S3")

#write.tree(all_tree$scenario.3, "sample.tre")

#Convert to pivot_longer() and filter only presences since ENFA
#is a presence-only model
spp_all_coords <- spp_all[, c(2, 3, 9:368)]
spp_all_coords <- spp_all_coords %>% pivot_longer(cols = 3:362, 
                                      values_to = "occurrence") %>% 
     filter(occurrence == 1)

#Adjusting names
spp_all_coords$name <- gsub("\\.", "_", spp_all_coords$name)

#Adjust names to filter only species that appeared on the tree
to_remove <- gsub(" ", "_", to_remove)
spp_all_coords2 <- spp_all_coords %>% filter(name %in% to_remove)
spp_all_coords2 <- as.data.frame(spp_all_coords2)

#Removing points that are still too close to each other,
#avoiding pseudoreplication (two presences on same cell)
spp_all_coords2 <- ecospat.occ.desaggregation2(xy = spp_all_coords2,
                                     min.dist = 0.125,
                                     by = "name")

#Remove species that after disaggregation got < 10 spatial records
#spp_all_coords_sum <- spp_all_coords2 %>% 
#     group_by(name) %>%
#     summarise(soma = sum(occurrence)) %>% 
#     filter(soma < 10) %>%
#     dplyr::select(name)
#to_filter <- as.data.frame((spp_all_coords_sum))
#to_filter$name <- as.character(to_filter$name)
#to_filter <- as.vector(to_filter$name)
#spp_all_coords <- spp_all_coords %>% filter(!name %in% to_filter)

#Prepare as spatial
spp_all_coords2 <- vect(spp_all_coords2, geom = c("x", "y"))

#Get EGV for each point
spp_all_coords2 <- terra::extract(preds, spp_all_coords2, bind = TRUE)
spp_all_coords2 <- terra::as.data.frame(spp_all_coords2, geom = "XY") #no NA

#Convert to sf
#spp_all_coords <- st_as_sf(x = spp_all_coords, 
#                  coords = c("x", "y"),
#                  crs = crs(preds))

all_names <- spp_all_coords2$name
all_names <- as.character(all_names)
spp_all_coords2 <- spp_all_coords2[, -1]

#Creating a list of species presence-background points
all_list <- split(x = spp_all_coords2, f = all_names)
#saveRDS(all_list, "rds/pre_all_list_cleaned_pivot_enfa.rds")
#all_list <- readRDS("rds/pre_all_list_cleaned_pivot_enfa.rds")

#SC Mask
sc <- vect("data/polygons/sc.shp")
new_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
sc <- project(sc, new_crs)

#Using all cells as background points
bg <- terra::xyFromCell(preds[[1]], 1:ncell(preds))
bg <- as.data.frame(bg)
bg <- vect(bg, geom = c("x", "y"))
bg <- terra::extract(preds, bg, bind = TRUE)
bg <- terra::as.data.frame(bg, geom = "XY")
bg <- na.omit(bg) #remove points that fall outside of the mask
bg$occurrence <- 0
bg <- bg[, c(8, 6, 7, 1:5)]

for(i in names(all_list)){
     all_list[[i]] <- rbind(all_list[[i]], bg)
                         }

#Adjusting input data to sf::data.frame format
for(i in names(all_list)){
     all_list[[i]] <- st_as_sf(x = all_list[[i]],
                               coords = c("x", "y"),
                               crs = crs(preds[[1]]))
}

#saveRDS(all_list, "rds/all_list_enfa.rds")

rm(list=ls())

################################################
################################################
########## MODELING PHYLOGENETIC ENFA ##########
################################################
################################################

#Loading species list
all_list <- readRDS("rds/all_list_enfa.rds")

#Loading predictors' extent/mask
preds <- readRDS("rds/predictors_df.rds")
preds <- rast(preds)
#crs(preds) <- NA
#The crs below also doesn't work
crs(preds) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#Loading phylogenetic tree
tree <- ape::read.tree("sample.tre")

#Modeling ENphylo
raremod <- ENphylo_modeling(input_data = all_list,
                            tree = tree,
                            input_mask = preds[[1]],
                            obs_col = "occurrence",
                            min_occ_enfa = 20,
                            boot_test_perc = 20,
                            boot_reps = 100,
                            eval_metric_for_imputation = "AUC",
                            output_options = "best",
                            clust = 1,
                            eval_threshold = 0.5)

#saveRDS(raremod, "rds/raremod.rds")

raremod <- readRDS("rds/raremod.rds")
main_dir <- getwd()

preds <- readRDS("rds/predictors_df.rds")
preds <- rast(preds)

#Predicting
raremod_pred <- ENphylo_prediction(object = raremod,
                                   newdata = preds,
                                   convert.to.suitability = T,
                                   output.dir = main_dir)


saveRDS(raremod_pred, "rds/raremod_pred.rds")

rm(list=ls())

#Reloading species modeled with and without ENFA
#Commons
spp_common <- terra::as.data.frame(readRDS("rds/pa_common_enfa.rds"),
                                   geom = "XY")
spp_common <- spp_common[, 2:279]
spp_common <- colnames(spp_common)
spp_common <- gsub("\\.", "_", spp_common)

#All (commons and rares)
enfapred <- readRDS("rds/raremod_pred.rds")
rares <- setdiff(names(enfapred), spp_common)

#Separating only rares (10-19) from ENphylo
enfa_rares <- list.files("ENphylo_prediction/",
                         full.names = TRUE,
                         pattern = paste0(rares, collapse = "|"))
enfa_rares <- rast(enfa_rares)
#Subsetting only suitability maps
enfa_rares <- enfa_rares["Suitability"]
enfa_rares <- as.list(enfa_rares)
names(enfa_rares) <- rares

#Crop on hfp extension/mask
hfp_zon <- rast("Zonation_data/condition/hfp2009.tif")

for(i in names(enfa_rares)){
     enfa_rares[[i]] <- terra::crop(x = enfa_rares[[i]], 
                                    y = hfp_zon, 
                                    mask = TRUE)
}


for(i in names(enfa_rares)){
     terra::writeRaster(rast(enfa_rares[i]), 
                        paste0("results/enfa_rares/", i, ".tif"),
                        overwrite = TRUE)
}

rares_stacked <- list.files("results/enfa_rares", full.names = TRUE, pattern = ".tif")
rares_stacked <- rast(rares_stacked)
rares_stacked <- round(sum(rares_stacked), 0)

#SSDM
writeRaster(rares_stacked, filename = "maps/enfa_rares_map.tif", overwrite = TRUE)



