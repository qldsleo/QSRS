cores <- as.numeric(args[1])
print(cores)

# load intersect soil points with covariates dataset
DSM_data <- read.csv(paste0(resd, '/Site_Covariate_intersect.csv'))
DSM_data <- DSM_data %>% drop_na(names(DSM_data)[10:ncol(DSM_data)])

# convert categorical variables to factors
for(fname_cov in covs_subset_cat){
  covname_Tmp <- gsub('.tif', '', fname_cov)
  if (covname_Tmp %in% names(DSM_data)){
    DSM_data[,covname_Tmp] <- as.factor(DSM_data[,covname_Tmp])
    print(paste0('Converted - ', covname_Tmp, ' - to factor'))
  }
}

# create directory for models to go in
dir.create(paste0(resd, method_dsm), showWarnings=F)
dir.create(paste0(resd, method_dsm, '/models'), showWarnings=F)
if(run_PCA){
  save_dir <- paste0(resd, method_dsm, '/models/PCA')
}else{
  save_dir <- paste0(resd, method_dsm, '/models/all')
}
dir.create(save_dir, showWarnings=F)

# initialise full df for saving results
names_valStats <- names(summaryStatsFn(0,0))
GOOFDat <- data.frame('Depths'=depths, stringsAsFactors=F) # Summary stats df
GOOFDat[,names_valStats] <- NA
# for validation data, also include PICP for the requested levels. 
GOOFDat[,paste0('PICP_', piclplot)] <- NA

### extra check of the reported 5th and 95th percentiles of prediction intervals
### the percentage where the validation data were less than the 5th percentile
### and the percentage where the validation data were greater than the 95th percentile.
### these are the two uncertainty indicators that get reported to give the prediciton interval - both stats should be around 5%
### if for eg data are very skewed, Q5 might often be too small and PLTQ5 might be > 5
### and Q95 might often be too large and PLTQ5 might be > 5
GOOFDat[,c('PLTQ5', 'PGTQ95')] <- NA

calGOOFDat <- data.frame('Depths'=depths, stringsAsFactors=F) # calibration data stats summary df
calGOOFDat[,paste0('Calib_', names_valStats)] <- NA

### data.frame for storing avgMSE for each depth layer
df_avgMSE <- data.frame('Depths'=depths, 'avgMSE'=NA, stringsAsFactors=F) 

### the quantiles of a predictive distribution to give the confidence levels given in piclplot - used in qrf method
q4pilb <- (1 - piclplot / 100) / 2
q4piub <- 1 - (1 - piclplot / 100) / 2

### the number of standard deviations to add (+/-) to the predicted value to give the required prediction intervals...
### used in the bootstrapping rf method
nsdsadd4q <- -1 * qnorm((1 - piclplot / 100) / 2)

# PCA feature reduction
# drops categorical variables first and adds them back on later
DSM_data.pca <- DSM_data %>% dplyr::select(-gsub('.tif', '', covs_subset_cat[which(gsub('.tif', '', covs_subset_cat) %in% names(DSM_data))]))
DSM_data.pca <- prcomp(na.omit(DSM_data.pca[,c(10:ncol(DSM_data.pca))]), center=T, scale.=T)
eigenvalues <- get_eigenvalue(DSM_data.pca)
write.csv(eigenvalues, file=paste0(resd, method_dsm, '/models/', '/PCA_eigenvalues.csv'))
dim <- eigenvalues[eigenvalues$cumulative.variance.percent > 95,] %>% rownames # might change to lower %
no_dim <- gsub('Dim.', '', dim[1]) %>% as.numeric
write.csv(DSM_data.pca$rotation, file=paste0(resd, method_dsm, '/models/', 'PCA_loadings.csv'))
loadings <- abs(DSM_data.pca$rotation)
loadings <- loadings[,1:no_dim]
pca_variables <- do.call(rbind, lapply(1:no_dim, function(i){
  rownames(loadings)[which.max(loadings[,i])]
})) %>% as.data.frame
pca_variables <- c(levels(as.factor(pca_variables$V1)), gsub('.tif','', covs_subset_cat[which(gsub('.tif', '', covs_subset_cat) %in% names(DSM_data))]))
write.csv(data.frame('covariates'=pca_variables), file=paste0(resd, method_dsm, '/models/', 'PCA_covariates.csv'))

## For each depth fit a model
# Establish calibration and validation data sets
print('Starting model fitting...')
cl <- makeCluster(cores)
registerDoSNOW(cl)
for (d in 1:length(depths)){
  print(depths[d])
  set.seed(123)
  
  if(identical(covs_subset, 'all')){
    covariates <- list.files(path=covd, pattern='.tif', full.names=F)
    covariates <- gsub('.tif', '', covariates[!covariates %in% covs_dropped])
  }else{
    covariates <- gsub('.tif', '', covs_subset)
  }
  
  if(d < 7){
    DSM_data_depth <- DSM_data[,c(1:3,3+d,10:ncol(DSM_data))]
    if(run_PCA){DSM_data_depth <- DSM_data_depth[,c(1:4,which(names(DSM_data_depth) %in% c(pca_variables)))]}
    DSM_data_depth <- DSM_data_depth[complete.cases(DSM_data_depth),]
    training_samples <- getTrainAndTestSamples_BalancedFactors(DSM_data_depth)
    training <- training_samples[[1]]
    cDat <- DSM_data_depth[training, 4:(ncol(DSM_data_depth))]
    if(!is.na(training_samples[[2]])){cDat <- cDat %>% dplyr::select(-c(training_samples[[2]]))}
    vDat <- DSM_data_depth[-training, 4:(ncol(DSM_data_depth))]
    if(!is.na(training_samples[[2]])){vDat <- vDat %>% dplyr::select(-c(training_samples[[2]]))}
  }else{
    DSM_data_depth <- DSM_data
    DSM_data_depth <- DSM_data_depth %>% 
      pivot_longer(-!c(X0.5.cm, X5.15.cm, X15.30.cm, X30.60.cm, X60.100.cm, X100.200.cm), names_to = 'Depth', values_to = 'VALUE') %>% as.data.frame
    DSM_data_depth <- DSM_data_depth %>% relocate(Depth, .after = ID) %>% relocate(VALUE, .after = ID)
    DSM_data_depth <- DSM_data_depth %>% mutate(Depth = recode(Depth, 'X0.5.cm'='2.5', 'X5.15.cm'='10', 'X15.30.cm'='22.5',
                                                               'X30.60.cm'='45', 'X60.100.cm'='80', 'X100.200.cm'='150'))
    DSM_data_depth$Depth <- as.numeric(DSM_data_depth$Depth)
    DSM_data_depth <- DSM_data_depth %>% drop_na
    if(run_PCA){DSM_data_depth <- DSM_data_depth[,c(1:5,which(names(DSM_data_depth) %in% c(pca_variables)))]}
    training_samples <- getTrainAndTestSamples_BalancedFactors(DSM_data_depth)
    training <- training_samples[[1]]
    x <- DSM_data_depth %>% dplyr::filter(Depth == 2.5)
    cDat <- DSM_data_depth[which(DSM_data_depth$ID %in% c(x$ID[training])),4:(ncol(DSM_data_depth))]
    vDat <- DSM_data_depth[which(DSM_data_depth$ID %in% c(x$ID[-training])),4:(ncol(DSM_data_depth))]
  }
  
  # uncertainty analysis using bootstrapping method
  dir.create(paste0(save_dir, '/', depths[d]), showWarnings=F) # Create models subdirectory for each depth
  
  fitControl <- trainControl(method='cv', number=k_kfold, returnResamp='final', verboseIter=F, indexFinal=seq(length(training))) #Ten K-folds

  TrainDat <- cDat
  names(TrainDat)[1] <- 'VALUE'
  
  if(method_dsm == 'qrf'){
    
    # run and save validation model
    fit_ranger <- train(VALUE ~ ., data=TrainDat, method='ranger', trControl=fitControl, importance='impurity', quantreg=T)
    modelFile <- paste0(save_dir, '/', depths[d], '/xvMod.rds')
    saveRDS(object=fit_ranger, file=modelFile)
    
    # fit model on full dataset for mapping
    fitControl_map <- trainControl(method='cv', number=k_kfold, returnResamp='final', verboseIter=F, indexFinal=seq(nrow(DSM_data_depth))) #Ten K-folds
    MapDat <- DSM_data_depth[,4:(ncol(DSM_data_depth))]
    names(MapDat)[1] <- 'VALUE'
    
    fit_ranger_map <- train(VALUE ~ ., data=MapDat, method='ranger', trControl=fitControl_map, importance='impurity', quantreg=T)
    modelFile_map <- paste0(save_dir, '/', depths[d], '/mapMod.rds')
    saveRDS(object=fit_ranger_map, file=modelFile_map)
    
    # calculate residuals
    xvResiduals <- data.frame('Y' = DSM_data_depth[training, 'Y'],
                            'X' = DSM_data_depth[training, 'X'],
                            'res' = resid(fit_ranger))
    write.csv(xvResiduals, paste0(save_dir, '/', depths[d], '/xvResiduals.csv'), row.names=F)
    
    mapResiduals <- data.frame('Y' = DSM_data_depth$Y,
                            'X' = DSM_data_depth$X,
                            'res' = resid(fit_ranger_map))
    write.csv(mapResiduals, paste0(save_dir, '/', depths[d], '/mapResiduals.csv'), row.names=F)

    # determine overall variable importance
    imp <- varImp(fit_ranger)
    impdf <- imp$importance
    if(d == 1){
      ovi <- impdf
      colnames(ovi)[d] <- depths[d]
    }else if (d == 2){
      ovi <- merge(ovi, impdf['Overall'], by='row.names', all.x=T)
      colnames(ovi)[d+1] <- depths[d]
    }else{
      impdf$Row.names <- row.names(impdf)
      ovi <- merge(ovi, impdf, by='Row.names', all.x=T)
      colnames(ovi)[d+1] <- depths[d]
    }
    
    # save all data, cDat and vDat data
    write.csv(cDat, paste0(save_dir, '/', depths[d], '/cDat.csv'), row.names=F)
    write.csv(vDat, paste0(save_dir, '/', depths[d], '/vDat.csv'), row.names=F)
    write.csv(MapDat, paste0(save_dir, '/', depths[d], '/mapDat.csv'), row.names=F)
    
    # assess goodness of fit
    # calibration data
    Dat <- cDat[1]
    Dat <- setNames(Dat, c('obs'))
    Dat$pred <- predict(fit_ranger, newdata=cDat)
    Dat$depth <- cDat$Depth
    
    if(d<7){
      calGOOFDat[d,paste0('Calib_', names_valStats)] <- summaryStatsFn(observed = Dat$obs, predicted = Dat$pred)
    }else{
      Dat_split <- split(Dat, Dat$depth)
      for(c in 1:length(Dat_split)){
        calGOOFDat[c+6,paste0('Calib_', names_valStats)] <- summaryStatsFn(observed = Dat_split[[c]]$obs, predicted = Dat_split[[c]]$pred)
        calGOOFDat$Depths[c+6] <- paste0('all_', calGOOFDat$Depths[c])
      }
    }

    if(d<7){cDat_GOOF <- data.frame('Depth' = depths[d], 'Obs' = Dat$obs, 'Pred' = Dat$pred)}else{
      cDat_GOOF <- data.frame('Depth' = cDat$Depth, 'Obs' = Dat$obs, 'Pred' = Dat$pred)
    }
    write.csv(cDat_GOOF, paste0(save_dir, '/', depths[d], '/cDat_GOOF.csv'), row.names=F)
    
    # validation data
    Dat <- vDat[1]
    Dat <- setNames( Dat, c('obs')) 
    Dat$pred <- predict(fit_ranger, newdata = vDat)
    Dat$depth <- vDat$Depth
    
    # pt_residuals <- st_as_sf(x = xvResiduals, coords = c("X", "Y"), crs = 'EPSG:4326')
    # pred_area <- st_as_sf(x = mapResiduals[-training,], coords = c("X", "Y"), crs = 'EPSG:4326') %>% as_Spatial %>% spTransform(CRS('+init=EPSG:3857'))
    # res_map <- autoKrige(res~1, input_data=as_Spatial(st_transform(pt_residuals, crs=proj4string(pred_area))), new_data=pred_area, model='Exp')
    # vDat_residuals <- res_map$krige_output[,"var1.pred"] %>% spTransform(CRS('+init=EPSG:4326')) %>% as.data.frame
    # vDat_residuals <- vDat_residuals$var1.pred
    # Dat$pred <- Dat$pred+vDat_residuals

    if(d<7){
      GOOFDat[d,names_valStats] <- summaryStatsFn(observed = Dat$obs, predicted = Dat$pred)
    }else{
      Dat_split <- split(Dat, Dat$depth)
      for(c in 1:length(Dat_split)){
        GOOFDat[c+6,names_valStats] <- summaryStatsFn(observed = Dat_split[[c]]$obs, predicted = Dat_split[[c]]$pred)
        GOOFDat$Depths[c+6] <- paste0('all_', GOOFDat$Depths[c])
      }
    }

    if(d<7){vDat_GOOF <- data.frame('Depth' = depths[d], 'Obs' = Dat$obs, 'Pred' = Dat$pred)}else{
      vDat_GOOF <- data.frame('Depth' = vDat$Depth, 'Obs' = Dat$obs, 'Pred' = Dat$pred)
    }
    write.csv(vDat_GOOF, paste0(save_dir, '/', depths[d], '/vDat_GOOF.csv'), row.names=F)
    
    # make each factor-level a column in vDat as an indicator
    for(fname_cov in covs_subset_cat){
      covname_Tmp <- gsub('.tif' , '' , fname_cov)
      if (covname_Tmp %in% names(fit_ranger$xlevels)){
        lvl_c <- fit_ranger$xlevels[[covname_Tmp]]
        for(lvl in lvl_c){
          vDat[,paste0(covname_Tmp, lvl)] <- as.numeric(vDat[,covname_Tmp] == lvl)
        }
      }
    }

    if (d<7){
      # lower and upper bounds of pred ints...
      lMat <- predict(fit_ranger$finalModel, vDat, type='quantiles', quantiles=q4pilb)$predictions
      uMat <- predict(fit_ranger$finalModel, vDat, type='quantiles', quantiles=q4piub)$predictions
      
      #PI coverage probability (PICP)
      # Determine which predictions are within the PI limits
      bMat <- matrix(NA, nrow=nrow(vDat), ncol=length(piclplot))
      for (i in 1:length(piclplot)) {
        names(vDat)[1] <- 'VALUE'
        bMat[, i] <- as.numeric(vDat$VALUE <= uMat[, i] & vDat$VALUE >= lMat[, i])
      }
      
      picpVec <- ((colSums(bMat)/nrow(bMat)) * 100)
      
      ### and the 5th 95th percentiles validation...
      i90 <- which(round(piclplot, digits = 1) == 90)
      if(length(i90) != 1){stop('Ouch - no piclplot defined for 90 percent prediction interval evaluation!')}
      
      b5 <- as.numeric(vDat$VALUE < lMat[, i90])
      pltq5 <- ((sum(b5)/length(b5)) * 100)
      
      b95 <- as.numeric(vDat$VALUE > uMat[, i90])
      pgtq95 <- ((sum(b95)/length(b95)) * 100)
      
      ### add coverage percentages to df of val stats...  
      GOOFDat[d,paste0('PICP_' , piclplot)] <- picpVec
      
      ### 5th and 95th percentiles validation...
      GOOFDat$PLTQ5[d] <- pltq5
      GOOFDat$PGTQ95[d] <- pgtq95
      
    }else{
      
      vDat_split <- split(vDat, vDat$Depth)
      
      for (x in 1:length(vDat_split)){
        # lower and upper bounds of pred ints...
        lMat <- predict(fit_ranger$finalModel, vDat_split[[x]], type='quantiles', quantiles=q4pilb)$predictions
        uMat <- predict(fit_ranger$finalModel, vDat_split[[x]], type='quantiles', quantiles=q4piub)$predictions
        
        #PI coverage probability (PICP)
        # Determine which predictions are within the PI limits
        bMat <- matrix(NA, nrow=nrow(vDat_split[[x]]), ncol=length(piclplot))
        for (i in 1:length(piclplot)) {
          names(vDat_split[[x]])[1] <- 'VALUE'
          bMat[, i] <- as.numeric(vDat_split[[x]]$VALUE <= uMat[, i] & vDat_split[[x]]$VALUE >= lMat[, i])
        }
        
        picpVec <- ((colSums(bMat)/nrow(bMat)) * 100)
        
        ### and the 5th 95th percentiles validation...
        i90 <- which(round(piclplot, digits = 1) == 90)
        if(length(i90) != 1){stop('Ouch - no piclplot defined for 90 percent prediction interval evaluation!')}
        
        b5 <- as.numeric(vDat_split[[x]]$VALUE < lMat[, i90])
        pltq5 <- ((sum(b5)/length(b5)) * 100)
        
        b95 <- as.numeric(vDat_split[[x]]$VALUE > uMat[, i90])
        pgtq95 <- ((sum(b95)/length(b95)) * 100)
        
        ### add coverage percentages to df of val stats...  
        GOOFDat[x+6,paste0('PICP_' , piclplot)] <- picpVec
        
        ### 5th and 95th percentiles validation...
        GOOFDat$PLTQ5[x+6] <- pltq5
        GOOFDat$PGTQ95[x+6] <- pgtq95
      }
    }
  } else if (method_dsm == 'mlm'){
    
    # run and save validation model
    for (i in 1:nbag){
      fit_mlm <- train(VALUE ~ ., data=TrainDat, method='lm', trControl=fitControl)
      modelFile <- paste0(resd, method_dsm, '/models/', ifelse(run_PCA,'PCA/','all/'), depths[d], '/xvMod_', i,'.rds')
      saveRDS(object = fit_mlm, file = modelFile)
    }
    
    r.models <- list.files(paste0(resd, method_dsm, '/models/', ifelse(run_PCA,'PCA/','all/'), depths[d]), pattern = "xvMod_", full.names=T)
    
    fit_mlm <- readRDS(r.models[1])
    imp <- varImp(fit_mlm)
    vi <- imp$importance
    for (i in 2:nbag) {
      fit_mlm <- readRDS(r.models[i])
      imp <- varImp(fit_mlm)
      impdf <- imp$importance
      vi <- cbind(vi, impdf)
    }
    vi$Mean <- rowMeans(vi) # Calculate mean across all bootstrap models
    
    mean <-  subset(vi, select = c('Mean'))
    
    if(d == 1){
      ovi <- mean
    }else{
      ovi <- cbind(ovi,mean)
    }
    colnames(ovi)[d] <- depths[d]

    # fit model on full dataset for mapping
    fitControl_map <- trainControl(method='cv', number=k_kfold, returnResamp='final', verboseIter=F, indexFinal=seq(nrow(DSM_data_depth))) #Ten K-folds
    MapDat <- DSM_data_depth[,4:(ncol(DSM_data_depth))]
    names(MapDat)[1] <- 'VALUE'
    
    fit_mlm_map <- train(VALUE ~ ., data=MapDat, method='lm', trControl=fitControl_map)
    modelFile_map <- paste0(resd, method_dsm, '/models/', ifelse(run_PCA,'PCA/','all/'), depths[d], '/mapMod.rds')
    saveRDS(object=fit_mlm_map, file=modelFile_map)
    
    # calculate residuals
    xvResiduals <- data.frame('Y' = DSM_data_depth[training, 'Y'],
                              'X' = DSM_data_depth[training, 'X'],
                              'res' = resid(fit_mlm))
    write.csv(xvResiduals, paste0(save_dir, '/', depths[d], '/xvResiduals.csv'), row.names=F)
    
    mapResiduals <- data.frame('Y' = DSM_data_depth$Y,
                               'X' = DSM_data_depth$X,
                               'res' = resid(fit_mlm_map))
    write.csv(mapResiduals, paste0(save_dir, '/', depths[d], '/mapResiduals.csv'), row.names=F)
    
    # save all data, cDat and vDat data
    write.csv(cDat, paste0(save_dir, '/', depths[d], '/cDat.csv'), row.names=F)
    write.csv(vDat, paste0(save_dir, '/', depths[d], '/vDat.csv'), row.names=F)
    write.csv(MapDat, paste0(save_dir, '/', depths[d], '/mapDat.csv'), row.names=F)
    
    # assess goodness of fit
    # calibration data
    Dat <- cDat[1]
    Dat <- setNames(Dat, c('obs'))
    Dat$pred <- predict(fit_mlm, newdata=cDat)
    Dat$depth <- cDat$Depth
    
    if(d<7){
      calGOOFDat[d,paste0('Calib_', names_valStats)] <- summaryStatsFn(observed = Dat$obs, predicted = Dat$pred)
    }else{
      Dat_split <- split(Dat, Dat$depth)
      for(c in 1:length(Dat_split)){
        calGOOFDat[c+6,paste0('Calib_', names_valStats)] <- summaryStatsFn(observed = Dat_split[[c]]$obs, predicted = Dat_split[[c]]$pred)
        calGOOFDat$Depths[c+6] <- paste0('all_', calGOOFDat$Depths[c])
      }
    }
    
    if(d<7){cDat_GOOF <- data.frame('Depth' = depths[d], 'Obs' = Dat$obs, 'Pred' = Dat$pred)}else{
      cDat_GOOF <- data.frame('Depth' = cDat$Depth, 'Obs' = Dat$obs, 'Pred' = Dat$pred)
    }
    write.csv(cDat_GOOF, paste0(save_dir, '/', depths[d], '/cDat_GOOF.csv'), row.names=F)
    
    # validation data
    Dat <- vDat[1]
    Dat <- setNames( Dat, c('obs')) 
    Dat$pred <- predict(fit_mlm, newdata = vDat)
    Dat$depth <- vDat$Depth
    
    # pt_residuals <- st_as_sf(x = xvResiduals, coords = c("X", "Y"), crs = 'EPSG:4326')
    # pred_area <- st_as_sf(x = mapResiduals[-training,], coords = c("X", "Y"), crs = 'EPSG:4326') %>% as_Spatial %>% spTransform(CRS('+init=EPSG:3857'))
    # res_map <- autoKrige(res~1, input_data=as_Spatial(st_transform(pt_residuals, crs=proj4string(pred_area))), new_data=pred_area, model='Exp')
    # vDat_residuals <- res_map$krige_output[,"var1.pred"] %>% spTransform(CRS('+init=EPSG:4326')) %>% as.data.frame
    # vDat_residuals <- vDat_residuals$var1.pred
    # Dat$pred <- Dat$pred+vDat_residuals
    
    if(d<7){
      GOOFDat[d,names_valStats] <- summaryStatsFn(observed = Dat$obs, predicted = Dat$pred)
    }else{
      Dat_split <- split(Dat, Dat$depth)
      for(c in 1:length(Dat_split)){
        GOOFDat[c+6,names_valStats] <- summaryStatsFn(observed = Dat_split[[c]]$obs, predicted = Dat_split[[c]]$pred)
        GOOFDat$Depths[c+6] <- paste0('all_', GOOFDat$Depths[c])
      }
    }
    
    if(d<7){vDat_GOOF <- data.frame('Depth' = depths[d], 'Obs' = Dat$obs, 'Pred' = Dat$pred)}else{
      vDat_GOOF <- data.frame('Depth' = vDat$Depth, 'Obs' = Dat$obs, 'Pred' = Dat$pred)
    }
    write.csv(vDat_GOOF, paste0(save_dir, '/', depths[d], '/vDat_GOOF.csv'), row.names=F)

    # make each factor-level a column in vDat as an indicator
    if(!is.null(covs_subset_cat)){
      for(fname_cov in covs_subset_cat){
        covname_Tmp <- gsub('.tif' , '' , fname_cov)
        if (covname_Tmp %in% names(fit_ranger$xlevels)){
          lvl_c <- fit_ranger$xlevels[[covname_Tmp]]
          for(lvl in lvl_c){
            vDat[,paste0(covname_Tmp, lvl)] <- as.numeric(vDat[,covname_Tmp] == lvl)
          }
        }else{}
      }
    }else{}
    
    if (d<7){
      # lower and upper bounds of pred ints...
      lMat <- predict(fit_ranger$finalModel, vDat, type='quantiles', quantiles=q4pilb)$predictions
      uMat <- predict(fit_ranger$finalModel, vDat, type='quantiles', quantiles=q4piub)$predictions

      #PI coverage probability (PICP)
      # Determine which predictions are within the PI limits
      bMat <- matrix(NA, nrow=nrow(vDat), ncol=length(piclplot))
      for (i in 1:length(piclplot)) {
        names(vDat)[1] <- 'VALUE'
        bMat[, i] <- as.numeric(vDat$VALUE <= uMat[, i] & vDat$VALUE >= lMat[, i])
      }
      
      picpVec <- ((colSums(bMat)/nrow(bMat)) * 100)
      
      ### and the 5th 95th percentiles validation...
      i90 <- which(round(piclplot, digits = 1) == 90)
      if(length(i90) != 1){stop('Ouch - no piclplot defined for 90 percent prediction interval evaluation!')}
      
      b5 <- as.numeric(vDat$VALUE < lMat[, i90])
      pltq5 <- ((sum(b5)/length(b5)) * 100)
      
      b95 <- as.numeric(vDat$VALUE > uMat[, i90])
      pgtq95 <- ((sum(b95)/length(b95)) * 100)
      
      ### add coverage percentages to df of val stats...  
      GOOFDat[d,paste0('PICP_' , piclplot)] <- picpVec
      
      ### 5th and 95th percentiles validation...
      GOOFDat$PLTQ5[d] <- pltq5
      GOOFDat$PGTQ95[d] <- pgtq95
      
    }else{
      
      vDat_split <- split(vDat, vDat$Depth)
      
      for (x in 1:length(vDat_split)){
        # lower and upper bounds of pred ints...
        lMat <- predict(fit_ranger$finalModel, vDat_split[[x]], type='quantiles', quantiles=q4pilb)$predictions
        uMat <- predict(fit_ranger$finalModel, vDat_split[[x]], type='quantiles', quantiles=q4piub)$predictions
        
        #PI coverage probability (PICP)
        # Determine which predictions are within the PI limits
        bMat <- matrix(NA, nrow=nrow(vDat_split[[x]]), ncol=length(piclplot))
        for (i in 1:length(piclplot)) {
          names(vDat_split[[x]])[1] <- 'VALUE'
          bMat[, i] <- as.numeric(vDat_split[[x]]$VALUE <= uMat[, i] & vDat_split[[x]]$VALUE >= lMat[, i])
        }
        
        picpVec <- ((colSums(bMat)/nrow(bMat)) * 100)
        
        ### and the 5th 95th percentiles validation...
        i90 <- which(round(piclplot, digits = 1) == 90)
        if(length(i90) != 1){stop('Ouch - no piclplot defined for 90 percent prediction interval evaluation!')}
        
        b5 <- as.numeric(vDat_split[[x]]$VALUE < lMat[, i90])
        pltq5 <- ((sum(b5)/length(b5)) * 100)
        
        b95 <- as.numeric(vDat_split[[x]]$VALUE > uMat[, i90])
        pgtq95 <- ((sum(b95)/length(b95)) * 100)
        
        ### add coverage percentages to df of val stats...  
        GOOFDat[x+6,paste0('PICP_' , piclplot)] <- picpVec
        
        ### 5th and 95th percentiles validation...
        GOOFDat$PLTQ5[x+6] <- pltq5
        GOOFDat$PGTQ95[x+6] <- pgtq95
      }
    }
    
  } else if (method_dsm == 'lmer'){
    ### lmer ###
    lmerft0 <- lmer(VALUE ~. , REML = FALSE , data = TrainDat) # for model comparison, use REML = FALSE
    lmerft <- step(lmerft0 , direction = 'back' , reduce.random = FALSE , alpha.fixed = 0.15) # alpha.fixed = 0.15 used to approx aic with 1 param
    lmerft <- get_model(lmerft)
    sumlmerft <- summary(lmerft)
  }
  
  else{
    stop('Error - unknown method_dsm!')
  }
}
stopCluster(cl)

# Summaries variable importance data
ovi$Mean <- rowMeans(ovi[2:7], na.rm=T) # Calculate mean across all depths

# ordering by mean importance.
ovi <- ovi[order(-ovi$Mean),] # Sort in descending order of mean importance

# Combine calibration and validation GOOF data
GOOFDat <- cbind(calGOOFDat, GOOFDat[,-1,drop=FALSE])

# Write Goodness of Fit data to file
write.csv(GOOFDat, paste0(save_dir, '/', 'GOOFData.csv'), row.names=F)

write.csv(ovi, paste0(save_dir, '/', 'Overall_variable_importance.csv'), row.names=F)

# plot Goodness of Fit
plot_data <- do.call(rbind, lapply(1:6, function(i){
  read.csv(paste0(save_dir, '/', depths[i], '/cDat_GOOF.csv'))
}))
plot_data$Depth = factor(plot_data$Depth, levels=c('X0.5.cm','X5.15.cm','X15.30.cm','X30.60.cm','X60.100.cm','X100.200.cm'))
ggplot(plot_data, aes(Obs, Pred))+
  geom_point(size=0.5)+geom_abline(col='red', size=0.75)+geom_smooth(method='lm',formula='y ~ x', col='blue', size=0.75)+
  facet_wrap(~ Depth, scales='free')+
  theme_bw()
ggsave(filename=paste0(save_dir, '/', 'Cal_alldepths.jpeg'))

plot_data <- do.call(rbind, lapply(1:6, function(i){
  read.csv(paste0(save_dir, '/', depths[i], '/vDat_GOOF.csv'))
}))
plot_data$Depth = factor(plot_data$Depth, levels=c('X0.5.cm','X5.15.cm','X15.30.cm','X30.60.cm','X60.100.cm','X100.200.cm'))
ggplot(plot_data, aes(Obs, Pred))+
  geom_point(size=0.5)+geom_abline(col='red', size=0.75)+geom_smooth(method='lm',formula='y ~ x', col='blue', size=0.75)+
  facet_wrap(~ Depth, scales='free')+
  theme_bw()
ggsave(filename=paste0(save_dir, '/','Val_alldepths.jpeg'))

if (length(depths) == 7){
  plot_data <- read.csv(paste0(save_dir, '/', depths[7], '/cDat_GOOF.csv'))
  plot_data$Depth = factor(plot_data$Depth, levels=c('2.5','10','22.5','45','80','150'))
  ggplot(plot_data, aes(Obs, Pred))+
    geom_point(size=0.5)+geom_abline(col='red', size=0.75)+geom_smooth(method='lm',formula='y ~ x', col='blue', size=0.75)+
    facet_wrap(~ Depth, scales='free', labeller=labeller(Depth=c('2.5'='X0.5.cm','10'='X5.15.cm','22.5'='X15.30.cm',
                                                                 '45'='X30.60.cm','80'='X60.100.cm','150'='X100.200.cm')))+
    theme_bw()
  ggsave(filename=paste0(save_dir, '/','all_Cal_alldepths.jpeg'))
  
  plot_data <- read.csv(paste0(save_dir, '/', depths[7], '/vDat_GOOF.csv'))
  plot_data$Depth = factor(plot_data$Depth, levels=c('2.5','10','22.5','45','80','150'))
  ggplot(plot_data, aes(Obs, Pred))+
    geom_point(size=0.5)+geom_abline(col='red', size=0.75)+geom_smooth(method='lm',formula='y ~ x', col='blue', size=0.75)+
    facet_wrap(~ Depth, scales='free', labeller=labeller(Depth=c('2.5'='X0.5.cm','10'='X5.15.cm','22.5'='X15.30.cm',
                                                                 '45'='X30.60.cm','80'='X60.100.cm','150'='X100.200.cm')))+
    theme_bw()
  ggsave(filename=paste0(save_dir, '/', 'all_Val_alldepths.jpeg'))
}

# PICP plots
GOOFDat_long <- GOOFDat[1:6,] %>% 
  pivot_longer(-!c(PICP_5, PICP_10, PICP_20, PICP_40, PICP_60, PICP_80, PICP_90, PICP_95, PICP_97.5, PICP_99), names_to='PICP', values_to='picpVec') %>% as.data.frame
GOOFDat_long$piclplot <- piclplot
GOOFDat_long$Depths = factor(GOOFDat_long$Depths, levels=c('X0.5.cm','X5.15.cm','X15.30.cm','X30.60.cm','X60.100.cm','X100.200.cm'))

ggplot(GOOFDat_long, aes(piclplot, picpVec))+
  geom_point(size=0.5)+geom_abline(col='red', size=0.75)+
  facet_wrap(~ Depths, scales='free')+
  theme_bw()+labs(y='PICP', x='Confidence level')
ggsave(filename=paste0(save_dir, '/', 'PICP_alldepths.jpeg'))

if (length(depths) == 7){
  GOOFDat_long <- GOOFDat[7:12,] %>% 
    pivot_longer(-!c(PICP_5, PICP_10, PICP_20, PICP_40, PICP_60, PICP_80, PICP_90, PICP_95, PICP_97.5, PICP_99), names_to='PICP', values_to='picpVec') %>% as.data.frame
  GOOFDat_long$piclplot <- piclplot
  GOOFDat_long$Depths <- gsub('all_', '', GOOFDat_long$Depths)
  GOOFDat_long$Depths = factor(GOOFDat_long$Depths, levels=c('X0.5.cm','X5.15.cm','X15.30.cm','X30.60.cm','X60.100.cm','X100.200.cm'))
  
  ggplot(GOOFDat_long, aes(piclplot, picpVec))+
    geom_point(size=0.5)+geom_abline(col='red', size=0.75)+
    facet_wrap(~ Depths, scales='free')+
    theme_bw()+labs(y='PICP', x='Confidence level')
  ggsave(filename=paste0(save_dir, '/', 'all_PICP_alldepths.jpeg'))
}

# create report
rmarkdown::render(input=rmarkdown,
                  output_file=paste0(attribute, '_v', version, '.html'),
                  output_dir=paste0(save_dir, '/'))