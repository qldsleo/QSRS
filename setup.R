###################################################################################################################################
###################################################################################################################################
###################################################      DSM set up     ###########################################################
###################################################################################################################################
###################################################################################################################################
# parsed arguments from batch scripts
args <- commandArgs(trailingOnly=T)

# set soil attribute and model version
attribute <- as.character(args[1])
version   <- as.numeric(args[2])
process   <- as.character(args[3]) # extract, model, model_sum, map, mosaic

# set directories and files
libpath         <- '/scratch/rsc3/leos/R/libs'
covd            <- '/scratch/rsc3/leos/QSRS/covariates'
extd            <- paste0('/scratch/rsc3/leos/QSRS/soil_attributes/', attribute)
resd            <- paste0('/scratch/rsc3/leos/QSRS/soil_attributes/', attribute, '/v', version, '/')
code            <- '/export/home/leos/QSRS/code/'
# fullfilename_td <- paste0(code, 'sali/', attribute, '/', attribute, '_harmonized.csv')
fullfilename_td <- '/scratch/rsc3/leos/QSRS/soil_attributes/all_SALI_data.csv'
fullfilename_A_inserts <- '/scratch/rsc3/leos/QSRS/soil_attributes/A_horizon_inserts.csv'
fullfilename_B_inserts <- '/scratch/rsc3/leos/QSRS/soil_attributes/B_horizon_inserts.csv'
slga_stats      <- '/scratch/rsc3/leos/QSRS/soil_attributes/slga_validation_stats.csv'

# modelling information
method_dsm     <- 'qrf' # qrf, cubist, mlm, rf
propn_training <- 0.7   # train:test split
k_kfold        <- 10    # number of k-folds
nbags          <- 10    # applicable to all models but qrf
piclplot       <- c(5, 10, 20, 40, 60, 80, 90, 95, 97.5, 99) # prediction interval confidence levels to validate
transform      <- 'NA' # 'NA' or 'log' or 'logit100' or 'sqrt'

# covariate versions
source(paste0(code, 'QSRS/versions.R'))

###################################################################################################################################
# load libraries
# install.packages("jpeg", lib=libpath)
.libPaths(libpath)
library(raster); library(terra); library(doParallel); library(ranger); library(caret); library(dplyr); library(tidyr); library(Metrics)
library(rmarkdown); library(factoextra); library(doSNOW); library(lubridate); library(mpspline2); library(DescTools); library(sf)
library(reshape); library(quickPlot); library(doMC); library(fastmap); library(plyr); library(sp); library(rgdal); library(Cubist)

# scripts
extraction <- paste0(code, 'QSRS/extraction.R')
modelling  <- paste0(code, 'QSRS/modelling.R')
model_sum  <- paste0(code, 'QSRS/modelling_summary.R')
plots      <- paste0(code, 'QSRS/plots.R')
mapping    <- paste0(code, 'QSRS/mapping.R')
mosaic     <- paste0(code, 'QSRS/mosaic.R')
rmarkdown  <- paste0(code, 'QSRS/QSRS.Rmd')
functions  <- paste0(code, 'QSRS/functions.R')
spline_clean <- paste0(code, 'QSRS/spline_clean.R')
source(functions)

# load training dataset
if (process == 'extract'){
  source(spline_clean)
} else {
  data <- read.csv(paste0(extd, '/harmonized.csv'))
}
print('Training data - ')
head(data)

# remove soil.depth column and check depths
if ('soil.depth' %in% colnames(data)){data <- data %>% dplyr::select(-soil.depth)}
depths <- c(names(data[-(1:3)])) # list of depths in training data
if(length(depths) > 6){depths <- depths[1:6]}
print(paste0('Depths - ', paste(depths, collapse=', ')))

# summarise data and check for outliers
summary(data)

###################################################################################################################################
# run selected process
if(process == 'extract'){
  source(extraction)
} else if (process == 'model'){
  source(modelling)
} else if(process == 'model_sum'){
  source(model_sum)
} else if (process == 'map'){
  source(mapping)
} else if (process == 'mosaic'){
  source(mosaic)
}