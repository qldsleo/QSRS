dep <- as.numeric(args[4])
print(dep)

# load intersect soil points with covariates dataset
DSM_data <- read.csv(paste0(extd, '/Site_Covariate_intersect.csv'))
if (!identical(covs_subset, 'all')){
  DSM_data <- DSM_data[,names(DSM_data) %in% c('Y', 'X', 'ID', depths, gsub('.tif', '', covs_subset))]
}
if (!identical(covs_dropped, '')){
  DSM_data <- DSM_data[,!names(DSM_data) %in% c(gsub('.tif', '', covs_dropped))]
}

# convert any categorical variables to factors
for(fname_cov in covs_subset_cat){
  covname_Tmp <- gsub('.tif', '', fname_cov)
  if (covname_Tmp %in% names(DSM_data)){
    DSM_data[,covname_Tmp] <- as.factor(DSM_data[,covname_Tmp])
    print(paste0('Converted - ', covname_Tmp, ' - to factor'))
  }
}

# create directory for models to go in
dir.create(resd, showWarnings=F)
dir.create(paste0(resd, method_dsm), showWarnings=F)
dir.create(paste0(resd, method_dsm, '/models'), showWarnings=F)
save_dir <- paste0(resd, method_dsm, '/models')

# model fitting
# cl <- makeCluster(cores); registerDoSNOW(cl)
for (d in dep){
  print(depths[d])
  set.seed(123)
  
  DSM_data_depth <- DSM_data[,c(1:3,3+d,(4+length(depths)):ncol(DSM_data))]
  DSM_data_depth <- DSM_data_depth[complete.cases(DSM_data_depth),]
  training_samples <- getTrainAndTestSamples_BalancedFactors(DSM_data_depth)
  training <- training_samples[[1]]
  cDat <- DSM_data_depth[training, 4:(ncol(DSM_data_depth))]
  if(!is.na(training_samples[[2]])){cDat <- cDat %>% dplyr::select(-c(training_samples[[2]]))}
  vDat <- DSM_data_depth[-training, 4:(ncol(DSM_data_depth))]
  if(!is.na(training_samples[[2]])){vDat <- vDat %>% dplyr::select(-c(training_samples[[2]]))}

  # create models subdirectory for each depth
  dir.create(paste0(save_dir, '/', depths[d]), showWarnings=F)
  
  # set up model training and tuning
  fitControl <- trainControl(method='cv', number=k_kfold, returnResamp='final', verboseIter=F, indexFinal=seq(length(training)), allowParallel=F)
  fitControl_map <- trainControl(method='cv', number=k_kfold, returnResamp='final', verboseIter=F, indexFinal=seq(nrow(DSM_data_depth)), allowParallel=F)
  
  TrainDat <- cDat
  names(TrainDat)[1] <- 'VALUE'
  MapDat <- DSM_data_depth[,4:(ncol(DSM_data_depth))]
  names(MapDat)[1] <- 'VALUE'
  
  # save all data, cDat and vDat data
  write.csv(cDat, paste0(save_dir, '/', depths[d], '/cDat.csv'), row.names=F)
  write.csv(vDat, paste0(save_dir, '/', depths[d], '/vDat.csv'), row.names=F)
  write.csv(MapDat, paste0(save_dir, '/', depths[d], '/mapDat.csv'), row.names=F)
  
  # transform data
  TrainDat$VALUE <- tfmFn(TrainDat$VALUE)
  MapDat$VALUE <- tfmFn(MapDat$VALUE)

########################################################################################################################################
  
  varimp_res_goof <- function(){
    
    # save variable importance
    r.models <- list.files(paste0(save_dir, '/', depths[d]), pattern = "xvMod_", full.names=T)
    fit_model <- readRDS(r.models[1])
    impdf <- varImp(fit_model)$importance
    for (i in 2:nbags){
      fit_model <- readRDS(r.models[i])
      vi <- varImp(fit_model)$importance
      impdf <- cbind(impdf, vi)
    }
    impdf$Mean <- rowMeans(impdf) # Calculate mean across all bootstrap models
    impdf <- subset(impdf, select = c('Mean'))
    colnames(impdf) <- depths[d]
    impdf$Covariate <- rownames(impdf)
    impdf <- impdf %>% relocate(Covariate)
    write.csv(impdf, paste0(save_dir, '/', depths[d], '/varImp.csv'), row.names=F)
    
    # assess goodness of fit - calibration and validation data
    r.models <- list.files(paste0(save_dir, '/', depths[d]), pattern = "xvMod_", full.names=T)
    fit_model <- readRDS(r.models[1])
    preddf <- predict(fit_model, newdata=cDat) %>% tfmFn(invt=T)
    for (i in 2:nbags){
      fit_model <- readRDS(r.models[i])
      pred <- predict(fit_model, newdata=cDat) %>% tfmFn(invt=T)
      preddf <- cbind(preddf, pred)
    }
    cDat_GOOF <- data.frame('Depth' = depths[d], 'Obs' = cDat[,1], 'Pred' = rowMeans(preddf))
    write.csv(cDat_GOOF, paste0(save_dir, '/', depths[d], '/cDat_GOOF.csv'), row.names=F)
    
    # need to double check this
    r.models <- list.files(paste0(save_dir, '/', depths[d]), pattern = "xvMod_", full.names=T)
    fit_model <- readRDS(r.models[1])
    preddf <- predict(fit_model, newdata=vDat) %>% tfmFn(invt=T)
    mseVec <- mse(vDat[,1], preddf)
    for (i in 2:nbags){
      fit_model <- readRDS(r.models[i])
      pred <- predict(fit_model, newdata=vDat) %>% tfmFn(invt=T)
      preddf <- cbind(preddf, pred)
      mseVec <- cbind(mseVec, mse(vDat[,1], pred))
    }
    vDat_GOOF <- data.frame('Depth' = depths[d], 'Obs' = vDat[,1], 'Pred' = rowMeans(preddf), 'Var' = apply(preddf, 1, var), 'MSE' = mean(mseVec))
    write.csv(vDat_GOOF, paste0(save_dir, '/', depths[d], '/vDat_GOOF.csv'), row.names=F)
    
  }
  
  if(method_dsm == 'qrf'){
    
    # fit model on training dataset for validation
    fit_model <- train(VALUE ~ ., data=TrainDat, method='ranger', trControl=fitControl, importance='impurity', quantreg=T)
    modelFile <- paste0(save_dir, '/', depths[d], '/xvMod.rds')
    saveRDS(object=fit_model, file=modelFile)
    
    # fit model on full dataset for mapping
    fit_model_map <- train(VALUE ~ ., data=MapDat, method='ranger', trControl=fitControl_map, importance='impurity', quantreg=T)
    modelFile_map <- paste0(save_dir, '/', depths[d], '/mapMod.rds')
    saveRDS(object=fit_model_map, file=modelFile_map)
    
    # save variable importance
    impdf <- varImp(fit_model)$importance
    colnames(impdf) <- depths[d]
    impdf$Covariate <- rownames(impdf)
    impdf <- impdf %>% relocate(Covariate)
    write.csv(impdf, paste0(save_dir, '/', depths[d], '/varImp.csv'), row.names=F)
    
    # lower and upper bounds of pred intervals
    for(fname_cov in covs_subset_cat){
      covname_Tmp <- gsub('.tif' , '' , fname_cov)
      if (covname_Tmp %in% names(fit_model$xlevels)){
        lvl_c <- fit_model$xlevels[[covname_Tmp]]
        for(lvl in lvl_c){
          vDat[,paste0(covname_Tmp, lvl)] <- as.numeric(vDat[,covname_Tmp] == lvl)
        }
      }
    }
    
    q4pilb <- (1 - piclplot / 100) / 2 # quantiles of a predictive distribution to give the confidence levels given in piclplot
    q4piub <- 1 - (1 - piclplot / 100) / 2 # quantiles of a predictive distribution to give the confidence levels given in piclplot
    
    lMat <- predict(fit_model$finalModel, vDat, type='quantiles', quantiles=q4pilb)$predictions %>% tfmFn(invt=T)
    write.csv(lMat, paste0(save_dir, '/', depths[d], '/lMat.csv'), row.names=F)
    uMat <- predict(fit_model$finalModel, vDat, type='quantiles', quantiles=q4piub)$predictions %>% tfmFn(invt=T)
    write.csv(uMat, paste0(save_dir, '/', depths[d], '/uMat.csv'), row.names=F)

    # assess goodness of fit - calibration and validation data
    cDat_GOOF <- data.frame('Depth' = depths[d], 'Obs' = cDat[,1], 'Pred' = predict(fit_model, newdata=cDat) %>% tfmFn(invt=T))
    write.csv(cDat_GOOF, paste0(save_dir, '/', depths[d], '/cDat_GOOF.csv'), row.names=F)

    vDat_GOOF <- data.frame('Depth' = depths[d], 'Obs' = vDat[,1], 'Pred' = predict(fit_model, newdata = vDat) %>% tfmFn(invt=T))
    write.csv(vDat_GOOF, paste0(save_dir, '/', depths[d], '/vDat_GOOF.csv'), row.names=F)
    
########################################################################################################################################
  } else if (method_dsm == 'mlm'){
    
    # fit model on training dataset for validation
    for (i in 1:nbags){
      fit_model <- train(VALUE ~ ., data=TrainDat, method='lm', trControl=fitControl)
      modelFile <- paste0(save_dir, '/', depths[d], '/xvMod_', i,'.rds')
      saveRDS(object = fit_model, file = modelFile)
    }
    
    # fit model on full dataset for mapping
    for (i in 1:nbags){
      fit_model_map <- train(VALUE ~ ., data=MapDat, method='lm', trControl=fitControl)
      modelFile_map <- paste0(save_dir, '/', depths[d], '/mapMod_', i,'.rds')
      saveRDS(object = fit_model_map, file = modelFile_map)
    }
    
    varimp_res_goof()
    
########################################################################################################################################
  } else if (method_dsm == 'cubist'){
    
    # fit model on training dataset for validation
    for (i in 1:nbags){
      fit_model <- train(VALUE ~ ., data=TrainDat, method='cubist', trControl=fitControl)
      modelFile <- paste0(save_dir, '/', depths[d], '/xvMod_', i,'.rds')
      saveRDS(object = fit_model, file = modelFile)
    }
    
    # fit model on full dataset for mapping
    for (i in 1:nbags){
      fit_model_map <- train(VALUE ~ ., data=MapDat, method='cubist', trControl=fitControl)
      modelFile_map <- paste0(save_dir, '/', depths[d], '/mapMod_', i,'.rds')
      saveRDS(object = fit_model_map, file = modelFile_map)
    }
    
    varimp_res_goof()
 
########################################################################################################################################
  } else if (method_dsm == 'rf'){
    
    # fit model on training dataset for validation
    for (i in 1:nbags){
      fit_model <- train(VALUE ~ ., data=TrainDat, method='ranger', trControl=fitControl, importance='impurity')
      modelFile <- paste0(save_dir, '/', depths[d], '/xvMod_', i,'.rds')
      saveRDS(object = fit_model, file = modelFile)
    }
    
    # fit model on full dataset for mapping
    for (i in 1:nbags){
      fit_model_map <- train(VALUE ~ ., data=MapDat, method='ranger', trControl=fitControl, importance='impurity')
      modelFile_map <- paste0(save_dir, '/', depths[d], '/mapMod_', i,'.rds')
      saveRDS(object = fit_model_map, file = modelFile_map)
    }
    
    varimp_res_goof()
    
########################################################################################################################################
  } else {
    stop('Error - unknown method_dsm!')
  }
}
# stopCluster(cl)
