chunk <- as.numeric(args[4])
print(chunk)
cores <- as.numeric(args[5])
print(cores)
dep <- as.numeric(args[6])
print(dep)
var <- as.numeric(args[7])
print(var)

# create directories
dir.create(paste0(resd, method_dsm, '/maps'), showWarnings=F)
save_map_dir <- paste0(resd, method_dsm, '/maps')

# stack covariates
if(identical(covs_subset, 'all')){
  list_ras <- list.files(path=covd, pattern='.tif', full.names=T)
  list_ras <- list_ras[!list_ras %in% paste0(covd, '/', covs_dropped)]
}else{
  list_ras <- paste0(covd, '/', covs_subset)
}
covariates.stack <- stack(list_ras)

# rename covariates in stack for consistency
if(identical(covs_subset, 'all')){
  list_ras <- list.files(path=covd, pattern='.tif', full.names=F)
  list_ras <- list_ras[!list_ras %in% covs_dropped]
}else{
  list_ras <- paste0(covs_subset)
}
names(covariates.stack) <- gsub('.tif','',list_ras)

# split into tiles on the fly
bs <- blockSize(covariates.stack, chunksize=2e10, minblocks=200)
i <- chunk
covariates.stack <- crop(covariates.stack, extent(covariates.stack, bs$row[[i]], bs$row[[i]]+bs$nrows[[i]], 1, ncol(covariates.stack)))

beginCluster(cores)
for (d in dep){
  
  # which depth
  pred_depth <- depths[d]
  
  # load model
  model <- readRDS(paste0(resd, method_dsm, '/models/', pred_depth, '/mapMod.rds'))
  
  # create directories
  dir.create(paste0(save_map_dir, '/', pred_depth), showWarnings=F)
  dir.create(paste0(save_map_dir, '/', pred_depth, '/pred'), showWarnings=F)
  dir.create(paste0(save_map_dir, '/', pred_depth, '/predvar05'), showWarnings=F)
  dir.create(paste0(save_map_dir, '/', pred_depth, '/predvar95'), showWarnings=F)
  
  # set filenames
  pred_file      <- paste0(save_map_dir, '/', pred_depth, '/pred/', i, '.tif')
  predvar05_file <- paste0(save_map_dir, '/', pred_depth, '/predvar05/', i, '.tif')
  predvar95_file <- paste0(save_map_dir, '/', pred_depth, '/predvar95/', i, '.tif')
  
  if (var == 50){
    clusterR(covariates.stack, predict, filename=pred_file, format="GTiff", overwrite=T, args=list(model=model))
  } else if (var == 5){
    t2 <- covariates.stack
    for(fname_cov in covs_subset_cat){
      covname_Tmp <- gsub('.tif', '', fname_cov)
      if (covname_Tmp %in% names(model$xlevels)){
        for(lvl in model$xlevels[[covname_Tmp]][-1]){
          t3 <- t2[[covname_Tmp]]
          t3[t3 != as.numeric(lvl)] <- 95157 # random number to ensure values are not changed in next line of code
          t3[t3 == as.numeric(lvl)] <- 1
          t3[t3 == 95157] <- 0
          names(t3) <- paste0(covname_Tmp, lvl)
          t2 <- addLayer(t2, t3)
        }
      }else{}
    }
  clusterR(t2, predict, filename=predvar05_file, format="GTiff", overwrite=T, args=list(model=model$finalModel, fun=function(model, ...) predict(model, type='quantiles', quantiles=c(0.05), ...)$predictions[,1]))
  } else if (var == 95){
    t2 <- covariates.stack
    for(fname_cov in covs_subset_cat){
      covname_Tmp <- gsub('.tif', '', fname_cov)
      if (covname_Tmp %in% names(model$xlevels)){
        for(lvl in model$xlevels[[covname_Tmp]][-1]){
          t3 <- t2[[covname_Tmp]]
          t3[t3 != as.numeric(lvl)] <- 95157 # random number to ensure values are not changed in next line of code
          t3[t3 == as.numeric(lvl)] <- 1
          t3[t3 == 95157] <- 0
          names(t3) <- paste0(covname_Tmp, lvl)
          t2 <- addLayer(t2, t3)
        }
      }else{}
    }
    clusterR(t2, predict, filename=predvar95_file, format="GTiff", overwrite=T, args=list(model=model$finalModel, fun=function(model, ...) predict(model, type='quantiles', quantiles=c(0.95), ...)$predictions[,1]))
  }
}
endCluster()
