##### Adjusting environmental data #####

#Unlink(".RData") if necessary
#remotes::install_github("rspatial/terra")


#Input filenames
inf <- list.files("data/rasters", pattern = "tif$", full.names = TRUE)

#Create output filenames and folder
outf <- gsub("data/rasters", "output", inf)
dir.create("output", FALSE, FALSE)

#Shapefile and extent to crop
sc <- vect("data/polygons/sc.shp")
new_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
sc <- project(sc, new_crs)
#sc <- buffer(sc, 10000)

#Aggregate to 5 km, crop, and mask
for (i in 1:length(inf)) {
     r <- terra::rast(inf[i])
     r <- terra::aggregate(r, 5, na.rm = TRUE)
     r <- terra::crop(r, sc, mask = TRUE, filename = outf[i], overwrite = TRUE)
                         }

##### Principal component analysis of environmental data #####
#Do the same for output/ and output_enfa/
vars <- list.files("output/", pattern = "tif$", full.name = TRUE)
vars <- rast(vars)

#Selecting only the variables
vars <- as.data.frame(vars, xy = TRUE)
vars <- na.omit(vars)

#Standardize variables
vars_scaled <- data.frame(apply(vars[, 3:25], 2, scale))

#PCA object: use return = TRUE to return the rotated values
pca_out <- prcomp(x = vars_scaled, retx = TRUE)

#Selecting components that represent ~95% of bioclimatic information
n_axes <- length(summary(pca_out)$importance[3, ])
cum_vars <- summary(pca_out)$importance[3, ]
vars_95 <- cum_vars <= 0.95 #five first axes

#Load scores
axes <- as.data.frame(pca_out$x) #rotated values
axes_95 <- axes[, vars_95] == T
axes_vars_95 <- axes[, 1:ncol(axes_95)]
axes_xy <- cbind(vars[, 1:2], axes_vars_95) #getting coordinates back
saveRDS(axes_xy, 'rds/predictors_df.rds') #Whole study region axes to be used

#As raster now
axes_xy_raster <- rast(axes_xy)
plot(axes_xy_raster)

#Adding a coordinate system
crs(axes_xy_raster) <- "epsg:4326"
saveRDS(axes_xy_raster, 'rds/predictors_raster.rds')
rm(list=ls())





