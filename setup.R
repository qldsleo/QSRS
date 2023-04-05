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
libpath         <- '/export/home/leos/libs'
covd            <- '/scratch/rsc3/leos/QSRS/covariates'
extd            <- paste0('/scratch/rsc3/leos/QSRS/soil_attributes/', attribute)
resd            <- paste0('/scratch/rsc3/leos/QSRS/soil_attributes/', attribute, '/v', version, '/')
code            <- '/export/home/leos/QSRS/code/'
fullfilename_td <- paste0(code, 'sali/', attribute, '/', attribute, '_harmonized.csv')

# modelling information
method_dsm     <- 'qrf' # qrf, cubist, mlm, rf
propn_training <- 0.7   # train:test split
k_kfold        <- 10    # number of k-folds
nbags          <- 10    # applicable to all models but qrf
piclplot       <- c(5, 10, 20, 40, 60, 80, 90, 95, 97.5, 99) # prediction interval confidence levels to validate

# covariates
covs_subset_cat <- c('Clim_GEOSS_Aus_Macroclimate_Bioclimatic_zones_qld.tif',
                     
                     'PM_GEOSS_Aus_Lithology_Lithology_qld.tif', 'PM_GEOSS_Aus_Lithology_Weathering_intensity', 'PM_Lithology_Map_Symbol_qld.tif',
                     'PM_Lithology_Min_Geol_Age_qld.tif', 'PM_Lithology_Unit_Type_qld.tif', 'PM_Silica_qld.tif',
                     
                     'Relief_Aus_Landform_Topographic_moisture_potential_qld.tif', 'Relief_GEOSS_Aus_Landform_Land_surface_forms_qld.tif',
                     'Relief_slope_relief.tif', 'Relief_tpi_class_1s_qld.tif', 'Relief_tpi_mask_1s_qld.tif',
                     
                     'Veg_GEOSS_Aus_Vegetation_structural_formations_qld.tif', 'Veg_IBRA_regions_qld.tif',
                     'Veg_LandCoverTrend_evi_class_qld.tif', 'Veg_preEuropeanVeg_qld.tif',
                     
                     'Landuse_qld.tif')

# include '.tif' if adding specific covariates
if (version == 1){
  covs_subset <- 'all'
  covs_dropped <- ''
  ###################################################################################################################################
} else if (version == 2){
  # dropped covariates which displayed paddock boundaries, roads etc as they were visible in final map (i.e. not realistic)
  covs_subset <- 'all'
  covs_dropped <- c('Clim_fwofs_qld.tif',
                    
                    'Other_BLUE_qld.tif', 'Other_CARBONATE_QUARTZ_BLUE_SWIR2_qld.tif',
                    'Other_HYDROXYL_1_PC2_qld.tif', 'Other_HYDROXYL_2_PC2_qld.tif', 'Other_ND_NIR_GREEN_qld.tif',
                    'Other_ND_RED_BLUE_qld.tif', 'Other_ND_RED_GREEN_qld.tif', 'Other_ND_SWIR1_BLUE_qld.tif',
                    'Other_ND_SWIR1_NIR_qld.tif', 'Other_ND_SWIR2_RED_qld.tif', 'Other_RED_qld.tif', 'Other_SWIR1_qld.tif',
                    
                    'Veg_FC_Max_BS_qld.tif', 'Veg_FC_Max_NPV_qld.tif',
                    'Veg_FC_Max_PV_qld.tif', 'Veg_FC_Mean_BS_qld.tif', 'Veg_FC_Mean_NPV_qld.tif',
                    'Veg_FC_Mean_PV_qld.tif', 'Veg_FC_Min_BS_qld.tif', 'Veg_FC_Min_NPV_qld.tif',
                    'Veg_FC_Min_PV_qld.tif', 'Veg_FC_SD_BS_qld.tif', 'Veg_FC_SD_NPV_qld.tif',
                    'Veg_FC_SD_PV_qld.tif', 'Veg_lztmre_aus_y20002011_dm7a2_d20050630_qld.tif', 'Veg_Persistant_green_Veg_qld.tif')
  ###################################################################################################################################
} else if (version == 3){
  # top 25 most important covariates from v1
  varimp <- read.csv(paste0(extd, '/v1/', method_dsm, '/models/Overall_variable_importance.csv'))
  covs_subset <- paste0(unique(gsub("\\d+$", "", varimp$Covariate[1:25])), '.tif')
  covs_dropped <- ''
  ###################################################################################################################################
} else if (version == 4){
  # top 25 most important covariates from v2
  varimp <- read.csv(paste0(extd, '/v2/', method_dsm, '/models/Overall_variable_importance.csv'))
  covs_subset <- paste0(unique(gsub("\\d+$", "", varimp$Covariate[1:25])), '.tif')
  covs_dropped <- ''
  ###################################################################################################################################
} else if (version == 5){
  # testing minimal covariates - consistent with most DSM papers + including land disturbance covariates for some attributes
  covs_subset <- c('Clim_PTA_qld.tif', 'Clim_TNI_qld.tif', 'Clim_TXX_qld.tif',
                   
                   'PM_Gravity_qld.tif', 'PM_radmap_v4_2019_filtered_pctk_GAPFilled_qld.tif',
                   'PM_radmap_v4_2019_filtered_ppmt_GAPFilled_qld.tif', 'PM_radmap_v4_2019_filtered_ppmu_GAPFilled_qld.tif',
                   'PM_Silica_qld.tif', 'PM_Weathering_Index_qld.tif',
                   
                   'Relief_twi_1s_qld.tif', 'Relief_slopepct1s_qld.tif', 'Relief_dem1sv1_0_qld.tif',
                   'Relief_mrrtf6g-a5_1s_qld.tif', 'Relief_mrvbf_int_qld.tif', 'Relief_plan_curvature_1s_qld.tif', 'Relief_PrescottIndex_01_1s_lzw_qld.tif',
                   
                   'Veg_preEuropeanVeg_qld.tif')
  covs_dropped <- ''
  if (any(attribute == c('tn','tp','ec','esp','phw','cat_ca','cat_cec','cat_k','cat_mg','cat_na','sar'))){
    covs_subset <- c(covs_subset, 'Veg_Persistant_green_Veg_qld.tif', 'Veg_LandCoverTrend_evi_mean_qld.tif')
  }
  ###################################################################################################################################
} else if (version == 6){
  # same as v5 but adding fractional cover
  covs_subset <- c('Clim_PTA_qld.tif', 'Clim_TNI_qld.tif', 'Clim_TXX_qld.tif',
                   
                   'PM_Gravity_qld.tif', 'PM_radmap_v4_2019_filtered_pctk_GAPFilled_qld.tif',
                   'PM_radmap_v4_2019_filtered_ppmt_GAPFilled_qld.tif', 'PM_radmap_v4_2019_filtered_ppmu_GAPFilled_qld.tif',
                   'PM_Silica_qld.tif', 'PM_Weathering_Index_qld.tif',
                   
                   'Relief_twi_1s_qld.tif', 'Relief_slopepct1s_qld.tif', 'Relief_dem1sv1_0_qld.tif',
                   'Relief_mrrtf6g-a5_1s_qld.tif', 'Relief_mrvbf_int_qld.tif', 'Relief_plan_curvature_1s_qld.tif', 'Relief_PrescottIndex_01_1s_lzw_qld.tif',
                   
                   'Veg_preEuropeanVeg_qld.tif', 'Veg_FC_Mean_BS_qld.tif', 'Veg_FC_Mean_NPV_qld.tif', 'Veg_FC_Mean_PV_qld.tif')
  covs_dropped <- ''
  if (any(attribute == c('tn','tp','ec','esp','phw','cat_ca','cat_cec','cat_k','cat_mg','cat_na','sar'))){
    covs_subset <- c(covs_subset, 'Veg_Persistant_green_Veg_qld.tif', 'Veg_LandCoverTrend_evi_mean_qld.tif')
  }
  ###################################################################################################################################
} else if (version == 7){}



###################################################################################################################################
# load libraries
# install.packages("jpeg", lib=libpath)
.libPaths(libpath)
library(raster); library(terra); library(doParallel); library(ranger); library(caret); library(dplyr); library(tidyr); library(Metrics)
library(rmarkdown); library(factoextra); library(doSNOW); library(lubridate); library(mpspline2); library(DescTools); library(sf)
library(reshape); library(quickPlot); library(doMC); library(fastmap); library(plyr); library(sp); library(rgdal); library(Cubist)

# load training dataset
data <- read.csv(fullfilename_td)
print('Training data - ')
head(data)

# remove soil.depth column and check depths
if ('soil.depth' %in% colnames(data)){data <- data %>% dplyr::select(-soil.depth)}
depths <- c(names(data[-(1:3)])) # list of depths in training data
if(length(depths) > 6){depths <- depths[1:6]} # this may need editing when doing 'all' depths
print(paste0('Depths - ', paste(depths, collapse=', ')))

# summarise data and check for outliers
summary(data)

# cleaned data
source(paste0(code, 'QSRS/clean.R'))
summary(data)

###################################################################################################################################
# scripts
extraction <- paste0(code, 'QSRS/extraction.R')
modelling  <- paste0(code, 'QSRS/modelling2.R')
model_sum  <- paste0(code, 'QSRS/modelling_summary.R')
plots      <- paste0(code, 'QSRS/plots.R')
mapping    <- paste0(code, 'QSRS/mapping.R')
mosaic     <- paste0(code, 'QSRS/mosaic.R')
rmarkdown  <- paste0(code, 'QSRS/QSRS.Rmd')
functions  <- paste0(code, 'QSRS/functions.R')

source(functions)

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