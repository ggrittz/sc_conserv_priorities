################################
################################
##### Adjusting field data #####
################################
################################

#Things are being done simultaneously for ENFA and SDM. Be careful!
#We did not proceed with ENFA analyses but the filenames were kept with ENFA names 
#for laziness sorry [22/04/2024].

#At the time of these analyses I wasn't so into GitHub so I didn't even use a version control,
#But everywhere is working as intended! [22/04/2024]

#Field data 1st survey coordinates for each site
coords <- read_xlsx('data/tables/UAs_IFFSC.xlsx', sheet = 1, skip = 1)
coords <- coords[, c("UA", "X (Longitude)", "Y (Latitude)")]
names(coords) <- c("ua", "x", "y")

#Field data survey site by species matrix
spp <- read_xlsx('data/tables/matriz_N_Ciclo1.xlsx', sheet = 1)

#Removing genus, subspecies, or dead ('morta')
spp$`Nome Científico` <- stringr::word(spp$`Nome Científico`, 1, 2)
spp <- na.omit(spp)

#More things to remove
to_remove <- c("Citrus X", "Cinnamomum sp.", "Mimosa sp.", "Ocotea sp.")
spp <- spp %>% filter(!`Nome Científico` %in% to_remove)

#Transpose on Excel and get it back
#write.csv(spp, 'data/tables/matrix_to_transpose.csv')

#Re-read the already transposed matrix
spp <- read.csv('data/tables/matrix_to_transpose.csv', sep = ';')

#Obtain a data set containing coordinates for each site plus site by species matrix
spp_coords <- merge(spp, coords, by = 'ua', all.x = TRUE, sort = FALSE)

#Convert abundance data to presence-absence in each sample unit
spp_coords[, 2:671] <- ifelse(spp_coords[, 2:671] > 0, 1, 0)

#Choosing species that have >= 20 individuals and also < 20
more_than_20 <- spp_coords[, 2:671] %>% dplyr::select(where(~sum(.) >= 20)) #278 species

less_than_20 <- spp_coords[, 2:671] %>% dplyr::select(where(~ sum(.) < 20)) #392 species

#Joining spatial info in each data set
more_than_20 <- cbind(spp_coords[, c("ua", "x", "y")], more_than_20)
less_than_20 <- cbind(spp_coords[, c("ua", "x", "y")], less_than_20)


############################
############################
##### Spatial cleaning #####
############################
############################

#Transforming each data set into spatial data
more_than_20 <- vect(more_than_20, geom = c("x", "y"), crs = "epsg:32722")
less_than_20 <- vect(less_than_20, geom = c("x", "y"), crs = "epsg:32722")

#Reproject both data sets to epsg:4326
more_than_20 <- project(x = more_than_20, y = "epsg:4326")
less_than_20 <- project(x = less_than_20, y = "epsg:4326")

#Getting back to data.frame format
common <- cbind(crds(more_than_20), as.data.frame(more_than_20@ptr$df$values()))
rare <- cbind(crds(less_than_20), as.data.frame(less_than_20@ptr$df$values()))

#Changing field names just for the sake
common$ua <- paste("ua_",common$ua, sep ="")
rare$ua <- paste("ua_",rare$ua, sep ="")

#Reading environmental data for extracting values at each survey site
#Do the same for enfa
env <- readRDS('rds/predictors_raster_enfa.rds')

#Transforming the data.frame to vector (points) to extract correctly
common <- vect(common, geom = c("x", "y"))
rare <- vect(rare, geom = c("x", "y"))

#Extracting and saving
common <- terra::extract(x = env, y = common, bind = TRUE, na.rm = TRUE)
crs(common) <- "epsg:4326"
#saveRDS(common, 'rds/pa_common_enfa.rds')

rare <- terra::extract(x = env, y = rare, bind = TRUE, na.rm = TRUE)
crs(rare) <- "epsg:4326"
#saveRDS(rare, 'rds/pa_rare_enfa.rds')


rm(list=ls())
gc()



