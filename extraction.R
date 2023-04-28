cores <- as.numeric(args[4])

# # check coord names
# data <- data %>% dplyr::rename(X = LONGITUDE, Y = LATITUDE)

# get covariate file names
files <- list.files(path=covd, pattern='.tif', full.names=T)

# stack covariate rasters
r1 <- stack(files)

# rename covariate rasters
names(r1) <- files %>% gsub(paste0(covd, '/'), '', .) %>% gsub('.tif', '', .)

# convert training data to points for extraction from intersecting covariate rasters
data_sp <- data
coordinates(data_sp) <- ~ X + Y
crs(data_sp) <- CRS('+init=EPSG:4326')
beginCluster(cores)
DSM_data <- raster::extract(r1, data_sp, sp=1, method='simple')
endCluster()
DSM_data <- as.data.frame(DSM_data)
DSM_data <- DSM_data %>% drop_na(names(DSM_data)[(3+length(depths)):ncol(DSM_data)])

# convert any categorical variables to factors
# covariates with factors NOT in QLD masked map are dropped entirely because we need them for mapping
for(fname_cov in covs_subset_cat){
  covname_Tmp <- gsub('.tif', '', fname_cov)
  if (!covname_Tmp %in% names(DSM_data)){
    print(paste0(fname_cov, ' - not in dataset'))
  }else{
    print(paste0('Checking - ', covname_Tmp, ' - factor levels'))
    DSM_data_levels <- as.numeric(levels(as.factor(na.omit(DSM_data[,covname_Tmp]))))
    print(paste0('Factors levels in training dataset - ', paste(DSM_data_levels, collapse=', ')))
    DSM_data[,covname_Tmp] <- as.factor(DSM_data[,covname_Tmp])
    covariate_raster <- rast(paste0(covd, '/', fname_cov)) # masked raster layer
    print(paste0('Min and max factor levels in - ', covname_Tmp, ' - ', paste0(minmax(covariate_raster)[,1], collapse=', ')))
    if (!all(minmax(covariate_raster)[,1] %in% DSM_data_levels)){
      DSM_data[,covname_Tmp] <- NULL
      print(paste0('WARNING! Dropping - ', covname_Tmp, ' - from dataset due to all factors in QLD rasters NOT in training dataset'))
    }else{
      covariate_levels <- terra::unique(covariate_raster)
      print(paste0('Factors levels in QLD raster - ', paste(covariate_levels[,1], collapse=', ')))
      if (!all(covariate_levels[,1] %in% DSM_data_levels)){
        DSM_data[,covname_Tmp] <- NULL
        print(paste0('WARNING! Dropping - ', covname_Tmp, ' - from dataset due to all factors in QLD rasters NOT in training dataset'))
      }
    }
  }
}

# create directory and save extraction
write.csv(DSM_data, file=paste0(extd, '/Site_Covariate_intersect.csv'), row.names=F)
