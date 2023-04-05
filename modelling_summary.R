# working directory
save_dir <- paste0(resd, method_dsm, '/models')

# initialise full df for saving results
names_valStats <- names(summaryStatsFn(0,0))
GOOFDat <- data.frame('Depths'=depths, stringsAsFactors=F) # validation data stats summary df
GOOFDat[,names_valStats] <- NA
GOOFDat[,paste0('PICP_', piclplot)] <- NA
GOOFDat[,c('PLTQ5', 'PGTQ95')] <- NA # % of validation data < 5th percentile & > 95th percentile - should be ~5%, otherwise data likely skewed
calGOOFDat <- data.frame('Depths'=depths, stringsAsFactors=F) # calibration data stats summary df
calGOOFDat[,paste0('Calib_', names_valStats)] <- NA
df_avgMSE <- data.frame('Depths'=depths, 'avgMSE'=NA, stringsAsFactors=F) # storing avgMSE for each depth layer
q4pilb <- (1 - piclplot / 100) / 2 # quantiles of a predictive distribution to give the confidence levels given in piclplot
q4piub <- 1 - (1 - piclplot / 100) / 2 # quantiles of a predictive distribution to give the confidence levels given in piclplot
nsdsadd4q <- -1 * qnorm((1 - piclplot / 100) / 2) # number of standard deviations to add (+/-) to the predicted value to give the required prediction intervals

for (d in 1:length(depths)){
  print(depths[d])

  cDat <- read.csv(paste0(save_dir, '/', depths[d], '/cDat.csv'))
  vDat <- read.csv(paste0(save_dir, '/', depths[d], '/vDat.csv'))
  
  if(method_dsm == 'qrf'){

    # determine overall variable importance
    ovi <- read.csv(paste0(save_dir, '/', depths[d], '/varImp.csv'))
    if (d > 1){impdf <- full_join(impdf, ovi)} else {impdf <- ovi}
    
    # assess goodness of fit - calibration and validation data
    Dat <- read.csv(paste0(save_dir, '/', depths[d], '/cDat_GOOF.csv'))
    calGOOFDat[d,paste0('Calib_', names_valStats)] <- summaryStatsFn(observed = Dat$Obs, predicted = Dat$Pred)

    Dat <- read.csv(paste0(save_dir, '/', depths[d], '/vDat_GOOF.csv'))
    GOOFDat[d,names_valStats] <- summaryStatsFn(observed = Dat$Obs, predicted = Dat$Pred)

    # PI coverage probability (PICP)
    # Determine which predictions are within the PI limits
    lMat <- read.csv(paste0(save_dir, '/', depths[d], '/lMat.csv'))
    uMat <- read.csv(paste0(save_dir, '/', depths[d], '/uMat.csv'))
    
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
      
########################################################################################################################################
  } else if (any(method_dsm == c('mlm','cubist','rf'))){
    
    # determine overall variable importance
    ovi <- read.csv(paste0(save_dir, '/', depths[d], '/varImp.csv'))
    if (d > 1){impdf <- full_join(impdf, ovi)} else {impdf <- ovi}
    
    # assess goodness of fit - calibration and validation data
    Dat <- read.csv(paste0(save_dir, '/', depths[d], '/cDat_GOOF.csv'))
    calGOOFDat[d,paste0('Calib_', names_valStats)] <- summaryStatsFn(observed = Dat$Obs, predicted = Dat$Pred)
    
    Dat <- read.csv(paste0(save_dir, '/', depths[d], '/vDat_GOOF.csv'))
    GOOFDat[d,names_valStats] <- summaryStatsFn(observed = Dat$Obs, predicted = Dat$Pred)

    # PI coverage probability (PICP)
    # Determine which predictions are within the PI limits
    val.sd <- sqrt(Dat$Var[1] + Dat$MSE[1])

    # zfactor multiplication
    vMat <- matrix(NA, nrow = nrow(vDat), ncol = length(piclplot))
    for (i in 1:length(piclplot)) {
      vMat[, i] <- val.sd * nsdsadd4q[i]
    }
    
    # Upper prediction limit
    uMat <- matrix(NA, nrow = nrow(vDat), ncol = length(piclplot))
    for (i in 1:length(piclplot)) {
      uMat[, i] <- Dat$Pred + vMat[, i]
    }
    
    # Lower prediction limit
    lMat <- matrix(NA, nrow = nrow(vDat), ncol = length(piclplot))
    for (i in 1:length(piclplot)) {
      lMat[, i] <- Dat$Pred - vMat[, i]
    }
    
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
    
  }
}

# Summaries variable importance data
impdf$Mean <- rowMeans(impdf[2:(1+length(depths))], na.rm=T)
ovi <- impdf[order(-impdf$Mean),]
write.csv(ovi, paste0(save_dir, '/', 'Overall_variable_importance.csv'), row.names=F)

# Combine calibration and validation GOOF data
GOOFDat <- cbind(calGOOFDat, GOOFDat[,-1,drop=FALSE])
write.csv(GOOFDat, paste0(save_dir, '/', 'GOOFData.csv'), row.names=F)

# plots
source(plots)

# example prediction
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

# crop stack
# qld <- st_read('/scratch/rsc3/leos/QSRS/QLD_regions/Local_Government_Areas.shp')
# bundaberg <- qld %>% dplyr::filter(ABBREV_NAM == 'BUNDABERG')
covariates.stack <- crop(covariates.stack, extent(152.1, 152.6, -25.08, -24.7))
model <- readRDS(paste0(resd, method_dsm, '/models/',depths[1],'/mapMod.rds'))
beginCluster(5)
prediction <- clusterR(covariates.stack, predict, args=list(model=model))
endCluster()

prediction2 <- crop(prediction, extent(152.295639, 152.358467, -24.981390, -24.921158))
prediction3 <- crop(prediction, extent(152.332203, 152.376835, -24.898894, -24.854510))

# create report
rmarkdown::render(input=rmarkdown,
                  output_file=paste0(attribute, '_v', version, '_', method_dsm, '.html'),
                  output_dir=paste0('/scratch/rsc3/leos/QSRS/soil_attributes/', attribute, '/'))
