library(terra)
#options(max.print = 10000)

sc_pa <- vect("data/polygons/sc_pas_cleaned.shp")
commons <- rast("Zonation_data/commons_out2/rankmap.tif")
rares <- rast("Zonation_data/rares_out2/rankmap.tif")

#Only values >= 0.7 to obtain 30% land conservation threshold

#Commons
commons_cbd <- clamp(commons, lower = 0.7, value = FALSE)
commons_cbd[commons_cbd > 0.7] <- 1
#How many cells?
length(cells(clamp(commons_cbd, lower = 0.7, value = FALSE))) #1613

#Rares
rares_cbd <- clamp(rares, lower = 0.7, value = FALSE)
rares_cbd[rares_cbd > 0.7] <- 1
#How many cells?
length(cells(clamp(rares_cbd, lower = 0.7, value = FALSE))) #1613, as expected (same res)

#What's the total nubmer of cells PAs can have?
#Using any raster here
cells(x = commons, y = sc_pa) #166 cells

#Union of commons and rares priority areas
union_cbd <- mosaic(commons_cbd, rares_cbd, fun = "max")
length(cells(union_cbd)) #2274 possible cells joining rares and commons priority areas

#Consensus (overlap) between commons and rares priority areas
intersect_cbd <- intersect(commons_cbd, rares_cbd)
intersect_cbd[intersect_cbd < 1] <- NA
table(values(intersect_cbd)) #952 cells by overlapping both maps

#How many of these cells are within PAs?
intersect_cbd_pa <- crop(intersect_cbd,
                         sc_pa,
                         mask = TRUE)
table(values(intersect_cbd_pa)) #94 consensus' cells are within PAs


#Are priority areas for common species more represented within PAs?
intersect_cbd_pa_commons <- crop(commons_cbd,
                              sc_pa,
                              mask = TRUE)
length(cells(intersect_cbd_pa_commons)) #142 common priority cells

#How about rares?
intersect_cbd_pa_rares <- crop(rares_cbd,
                               sc_pa,
                               mask = TRUE)
length(cells(intersect_cbd_pa_rares)) #114 cells

#Saving all products for plotting in QGIS
writeRaster(commons_cbd, "Zonation_data/Results/commons_cbd.tif")
writeRaster(rares_cbd, "Zonation_data/Results/rares_cbd.tif")
writeRaster(intersect_cbd, "Zonation_data/Results/intersect_commons_rares_cbd.tif")
writeRaster(union_cbd, "Zonation_data/Results/union_commons_rares_cbd.tif")












