dep <- as.numeric(args[1])
print(dep)

# which depth
d <- dep
pred_depth <- depths[d]

map_dir <- paste0(resd, method_dsm, '/maps/', ifelse(run_PCA,'PCA/','all/'))
no_files <- list.files(paste0(map_dir, '/', pred_depth, '/pred'), full.names=T)

if(length(no_files) != 200){stop('There are not 200 files for mosaicing. Check mapping process')}
#####################################
##### MERGE MEAN PREDICTED MAPS #####
rast.list <- list()
for (i in 1:length(no_files)){
  tile_temp <- raster(paste0(map_dir, '/', pred_depth, '/pred/', i, '.tif'))
  rast.list[[i]] <- tile_temp
}

names(rast.list)[1:2] <- c('x', 'y')
rast.list$fun <- mean
rast.list$na.rm <- TRUE
rast.list$tolerance <- 0.1
y_pred <- do.call(raster::mosaic, rast.list)

depth <- gsub('X', '', pred_depth)
depth <- gsub('.cm', '', depth)
depth <- gsub('[.]', '-', depth)
writeRaster(y_pred, filename=paste0(map_dir, attribute, '_mean_', depth, 'cm_30m_v', version, '.tif'), overwrite=T)

# Plot(y_pred, title=paste0(attribute, ' - ', pred_depth))

#####################################
### MERGE LOWER CI PREDICTED MAPS ###
# rast.list <- list()
# for (i in 1:length(no_files)){
#   tile_temp <- raster(paste0(resd, 'maps/', ifelse(run_PCA,'PCA/','all/'), pred_depth, '/predvar05/', i, '.tif'))
#   rast.list[[i]] <- tile_temp
# }
# 
# names(rast.list)[1:2] <- c('x', 'y')
# rast.list$fun <- mean
# rast.list$na.rm <- TRUE
# rast.list$tolerance <- 0.1
# y_pred <- do.call(raster::mosaic, rast.list)
# 
# depth <- gsub('X', '', pred_depth)
# writeRaster(y_pred, filename=paste0(resd, 'maps/', ifelse(run_PCA,'PCA/','all/'), attribute, '_var05_', depth, 'cm_30m_v', version, '.tif'), overwrite=T)
# 
# #####################################
# ### MERGE UPPER CI PREDICTED MAPS ###
# rast.list <- list()
# for (i in 1:length(no_files)){
#   tile_temp <- raster(paste0(resd, 'maps/', ifelse(run_PCA,'PCA/','all/'), pred_depth, '/predvar95/', i, '.tif'))
#   rast.list[[i]] <- tile_temp
# }
# 
# names(rast.list)[1:2] <- c('x', 'y')
# rast.list$fun <- mean
# rast.list$na.rm <- TRUE
# rast.list$tolerance <- 0.1
# y_pred <- do.call(raster::mosaic, rast.list)
# 
# depth <- gsub('X', '', pred_depth)
# writeRaster(y_pred, filename=paste0(resd, 'maps/', ifelse(run_PCA,'PCA/','all/'), attribute, '_var95_', depth, 'cm_30m_v', version, '.tif'), overwrite=T)
