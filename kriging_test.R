libpath <- '/export/home/leos/libs'; .libPaths(libpath)
library(sp)
library(gstat)
library(automap)
library(raster)
library(sf)

d <- 1

residuals <- read.csv(paste0(resd, '/models/', ifelse(run_PCA,'PCA/','all/'), depths[d], '/residuals.csv'))
residuals <- na.omit(residuals)
residuals <- st_as_sf(x = residuals, coords = c("X", "Y"), crs = 'EPSG:4326')

lzn.fit <- autofitVariogram(res~1, input_data=as_Spatial(residuals), model=c('Exp'))
plot(lzn.fit)

pred_area <- covariates.stack[[1]] %>% raster %>% as('SpatialPoints') %>% spTransform(CRS('+init=EPSG:3857'))

res_map <- autoKrige(res~1, input_data=as_Spatial(st_transform(residuals, crs=proj4string(pred_area))), new_data=pred_area, model='Sph')

krig_map <- rasterFromXYZ(res_map$krige_output[,"var1.pred"])

