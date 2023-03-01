###################################################################################################################################
###################################################################################################################################
###################################################      DSM set up     ###########################################################
###################################################################################################################################
###################################################################################################################################

# set library location
libpath <- '/export/home/leos/libs'

# set soil attribute and model version
# clay, silt, fs, cs, cat_ca, cat_mg, cat_na, cat_k, cat_cec, tn, tp (others: clay_act, esp, sar)
# need to extract ---- ec, phw, col_p, wb_oc
attribute <- 'clay'
version   <- 3

# reduce number covariates with a PCA?
run_PCA       <- FALSE    # will run a PCA to reduce covariates - if using PCA, leave TRUE throughout entire process

# what to process?
extract_now   <- FALSE   # step 1: will extract point-covariate intersecting data
xv_now        <- FALSE   # step 2: will run modelling script
map_now       <- TRUE   # step 3: will run mapping script
mosaic_now    <- FALSE   # step 4: will run mosaic script
cleanup_now   <- FALSE   # step 5: will remove 200 tiles generated from mapping script - only use after mosaicing is completed

# set directories and files
covd            <- '/scratch/rsc3/leos/QSRS/covariates'
resd            <- paste0('/scratch/rsc3/leos/QSRS/soil_attributes/', attribute, '/v', version, '/')
code            <- '/export/home/leos/QSRS/code/'
fullfilename_td <- paste0(code, 'sali/', attribute, '/', attribute, '_harmonized.csv')

# modelling information
method_dsm     <- 'qrf' # qrf, cubist, mlm
propn_training <- 0.7   # train:test split
k_kfold        <- 10    # number of k-folds
nbags          <- 50    # applicable to all models but qrf
piclplot       <- c(5, 10, 20, 40, 60, 80, 90, 95, 97.5, 99) # prediction interval confidence levels to validate

# covariates
covs_subset_cat <- c('Clim_GEOSS_Aus_Macroclimate_Bioclimatic_zones_qld.tif', 'PM_GEOSS_Aus_Lithology_Lithology_qld.tif',
                     'PM_Lithology_Map_Symbol_qld.tif', 'PM_Lithology_Unit_Type_qld.tif',
                     'Relief_Aus_Landform_Topographic_moisture_potential_qld.tif', 'Relief_GEOSS_Aus_Landform_Land_surface_forms_qld.tif',
                     'Relief_tpi_class_1s_qld.tif', 'Relief_tpi_mask_1s_qld.tif',
                     'Veg_GEOSS_Aus_Vegetation_structural_formations_qld.tif', 'Veg_IBRA_regions_qld.tif',
                     'Veg_LandCoverTrend_evi_class_qld.tif', 'Veg_preEuropeanVeg_qld.tif',
                     'Soil_QSRS_qld.tif')

if (version == 1){
  covs_subset <- 'all' # if subsetting, include '.tif'
  covs_dropped <- ''
} else if (version == 2){
  covs_subset <- 'all' # if subsetting, include '.tif'
  covs_dropped <- c('Clim_fwofs_qld.tif', 'Other_BLUE_qld.tif', 'Other_CARBONATE_QUARTZ_BLUE_SWIR2_qld.tif',
                    'Other_HYDROXYL_1_PC2_qld.tif', 'Other_HYDROXYL_2_PC2_qld.tif', 'Other_ND_NIR_GREEN_qld.tif',
                    'Other_ND_RED_BLUE_qld.tif', 'Other_ND_RED_GREEN_qld.tif', 'Other_ND_SWIR1_BLUE_qld.tif',
                    'Other_ND_SWIR1_NIR_qld.tif', 'Other_ND_SWIR2_RED_qld.tif', 'Other_RED_qld.tif',
                    'Other_SWIR1_qld.tif', 'Veg_FC_Max_BS_qld.tif', 'Veg_FC_Max_NPV_qld.tif',
                    'Veg_FC_Max_PV_qld.tif', 'Veg_FC_Mean_BS_qld.tif', 'Veg_FC_Mean_NPV_qld.tif',
                    'Veg_FC_Mean_PV_qld.tif', 'Veg_FC_Min_BS_qld.tif', 'Veg_FC_Min_NPV_qld.tif',
                    'Veg_FC_Min_PV_qld.tif', 'Veg_FC_SD_BS_qld.tif', 'Veg_FC_SD_NPV_qld.tif',
                    'Veg_FC_SD_PV_qld.tif', 'Veg_lztmre_aus_y20002011_dm7a2_d20050630_qld.tif', 'Veg_Persistant_green_Veg_qld.tif')
} else if (version == 3){
  covs_subset_cat <- '' # 'Soil_QSRS_qld.tif'
  covs_subset <- c('Clim_PTA_qld.tif', 'Clim_TNI_qld.tif', 'Clim_TXX_qld.tif', # climate
                   'PM_radmap_v4_2019_filtered_pctk_GAPFilled_qld.tif', 'PM_radmap_v4_2019_filtered_ppmt_GAPFilled_qld.tif', # parent material
                   'PM_radmap_v4_2019_filtered_ppmu_GAPFilled_qld.tif', 'PM_Silica_qld.tif',
                   'PM_Weathering_Index_qld.tif', # parent material
                   # 'Soil_QSRS_qld.tif', # soil
                   'Relief_twi_1s_qld.tif', 'Relief_slopepct1s_qld.tif', 'Relief_dem1sv1_0_qld.tif', # relief
                   'Relief_mrvbf_int_qld.tif', 'Relief_mrrtf6g-a5_1s_qld.tif', 'Relief_PrescottIndex_01_1s_lzw_qld.tif' # relief
                   # 'Veg_FC_Mean_BS_qld.tif', 'Veg_FC_Mean_PV_qld.tif', 'Veg_FC_Mean_NPV_qld.tif' # organisms
                   )
  covs_dropped <- NA
}

###################################################################################################################################
# load libraries
# install.packages("qrnn", lib=libpath)
.libPaths(libpath)
library(raster); library(terra); library(doParallel); library(ranger); library(caret); library(dplyr); library(tidyr)
library(rmarkdown); library(factoextra); library(doSNOW); library(lubridate); library(mpspline2); library(DescTools)
library(reshape); library(quickPlot); library(doMC); library(fastmap); library(plyr); library(sp); library(rgdal); library(Cubist)

# load training dataset
data <- read.csv(fullfilename_td)

# check coord names
if((toupper(names(data)[1]) == 'LONGITUDE') & (toupper(names(data)[2]) == 'LATITUDE')){
  names(data)[1] <- 'X'
  names(data)[2] <- 'Y'
}else if((toupper(names(data)[2]) == 'LONGITUDE') & (toupper(names(data)[1]) == 'LATITUDE')){
  names(data)[1] <- 'Y'
  names(data)[2] <- 'X'
}else{
  stop('Stopping - not sure which are the X/Y coordinates in data file!')
}

print('Training data - ')
head(data)

# remove soil.depth column and check depths
if ('soil.depth' %in% colnames(data)){data <- data %>% dplyr::select(-soil.depth)}
depths <- c(names(data[-(1:3)])) # list of depths in training data
# if(length(depths) > 6){depths <- depths[1:6]}else{} # this may need editing when doing 'all' depths
print(paste0('Depths - ', paste(depths, collapse=', ')))

# summarise data and check for outliers
print('Training data without cleaning - ')
summary(data)

# clean dataset
for (i in 4:9){
  if (attribute == 'clay' | attribute == 'cs' | attribute == 'fs' | attribute == 'silt'){
    data[i] <- replace(data[i], data[i] < 0, NA)
    data[i] <- replace(data[i], data[i] > 100, NA)
  } else if (attribute == 'phw' | attribute == 'phcl'){
    data[i] <- replace(data[i], data[i] < 0, NA)
    data[i] <- replace(data[i], data[i] > 14, NA)
  } else if (attribute == 'cat_ca'){
    data[i] <- replace(data[i], data[i] > 200, NA)
  } else if (attribute == 'cat_cec'){
    data[i] <- replace(data[i], data[i] > 300, NA)
  } else if (attribute == 'cat_k'){
    data[i] <- replace(data[i], data[i] > 15, NA)
  } else if (attribute == 'cat_na'){
    data[i] <- replace(data[i], data[i] > 150, NA)
  } else if (attribute == 'cat_mg'){
    data[i] <- replace(data[i], data[i] > 100, NA)
  } else if (attribute == 'tn'){
    data[i] <- replace(data[i], data[i] > 3, NA)
  } else if (attribute == 'tp'){
    data[i] <- replace(data[i], data[i] > 2, NA)
  }
}

print('Training data with cleaning - ')
summary(data)

###################################################################################################################################
# scripts
extraction <- paste0(code, 'QSRS/extraction.R')
modelling  <- paste0(code, 'QSRS/modelling.R')
mapping    <- paste0(code, 'QSRS/mapping.R')
mosaic     <- paste0(code, 'QSRS/mosaic.R')
rmarkdown  <- paste0(code, 'QSRS/QSRS.Rmd')
functions  <- paste0(code, 'QSRS/functions.R')

# parsed arguments from batch scripts
args <- commandArgs(trailingOnly=T)

source(functions)

if(extract_now){source(extraction)}

if(xv_now){source(modelling)}

if(map_now){source(mapping)}

if(mosaic_now){source(mosaic)}