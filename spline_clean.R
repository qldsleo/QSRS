# read in csv
extract <- read.csv(fullfilename_td)
extract <- extract %>% mutate(sand = cs+fs,
                              sand_method = paste0(cs_method,'_',fs_method))
extract <- extract %>% select(ID, UD, LD, paste0(attribute), paste0(attribute, '_method'), OBJECTID, PROJECT_CODE, SITE_ID, OBS_NO, HORIZON_NO, SAMPLE_NO,
                              YEAR, LATITUDE, LONGITUDE, ZONE, EASTING, NORTHING, DISTURB_TYPE, NATURE)
names(extract)[which(names(extract) == paste0(attribute))] <- "VALUE"
extract <- extract %>% filter(!is.na(VALUE))
extract <- distinct(extract, ID, UD, LD, .keep_all=T)

## copied from https://github.com/bitproxima/Modelling-mapping-with-Cubist/blob/master/Step1_SiteData.sql#L292
# standard filters
extract <- extract %>% filter(!is.na(VALUE),
                                  VALUE > 0,
                                  !(PROJECT_CODE == 'CQC'),
                                  !(PROJECT_CODE == 'EIM' & SITE_ID == 6051),
                                  !(PROJECT_CODE == 'QCS' & SITE_ID %in% c(20, 21, 85, 86)),
                                  !(PROJECT_CODE == 'FSE' & SITE_ID == 126 & SAMPLE_NO == 3),
                                  !(PROJECT_CODE == 'ABC' & SITE_ID == 315))

# specific filters for each attribute
if (attribute == 'clay' | attribute == 'cs' | attribute == 'fs' | attribute == 'sand'){
  extract <- extract %>% filter(VALUE <= 100,
                                !(PROJECT_CODE == 'BAN' & SITE_ID == 95),
                                !(PROJECT_CODE == 'BAMAR' & SITE_ID == 952 & SAMPLE_NO == 5))
} else if (attribute == 'esp'){
  extract <- extract %>% filter(VALUE <= 100,
                                !(PROJECT_CODE == 'ABC' & SITE_ID == 505  & SAMPLE_NO == 33),
                                !(PROJECT_CODE == 'ABC' & SITE_ID == 505  & SAMPLE_NO == 36),
                                !(PROJECT_CODE == 'ABC' & SITE_ID == 500  & SAMPLE_NO == 31),
                                !(PROJECT_CODE == 'WDH' & SITE_ID == 9098 & SAMPLE_NO == 13))
} else if (attribute == 'ec'){
  extract <- extract %>% filter(!(PROJECT_CODE == 'MON'   & SITE_ID == 6094 & SAMPLE_NO == 4),
                                !(PROJECT_CODE == 'MON'   & SITE_ID == 6089 & SAMPLE_NO == 3),
                                !(PROJECT_CODE == 'AGOPS' & SITE_ID == 162  & SAMPLE_NO == 16),
                                !(PROJECT_CODE == 'EIL'   & SITE_ID == 1000))
} else if (attribute == 'silt'){
  extract <- extract %>% filter(VALUE <= 100,
                                !(PROJECT_CODE == 'BAN'   & SITE_ID == 95),
                                !(PROJECT_CODE == 'MCL'   & SITE_ID == 9052 & SAMPLE_NO == 31),
                                !(PROJECT_CODE == 'BAMAR' & SITE_ID == 952  & SAMPLE_NO == 5),
                                !(PROJECT_CODE == 'MON'),
                                !(PROJECT_CODE == 'CCL'   & SITE_ID == 317  & SAMPLE_NO == 2))
} else if (attribute == 'cat_ca'){
  extract <- extract %>% filter(!(PROJECT_CODE == 'BAMAR' & SITE_ID == 952  & SAMPLE_NO == 5),
                                !(PROJECT_CODE == 'EIR'   & SITE_ID == 9021 & SAMPLE_NO == 36),
                                !(PROJECT_CODE == 'SALTC' & SITE_ID == 400  & SAMPLE_NO == 3),
                                !(PROJECT_CODE == 'SALTC' & SITE_ID == 400  & SAMPLE_NO == 4),
                                !(PROJECT_CODE == 'BDSM'  & SITE_ID == 345  & SAMPLE_NO == 1),
                                !(PROJECT_CODE == '3MC'   & SITE_ID == 9014))
} else if (attribute == 'cat_mg'){
  extract <- extract %>% filter(!(PROJECT_CODE == 'BDSM' & SITE_ID == 345 & SAMPLE_NO == 1))
} else if (attribute == 'cat_na'){
  extract <- extract %>% filter(!(PROJECT_CODE == 'CQA' & SITE_ID == 1001),
                                !(PROJECT_CODE == 'ABC' & SITE_ID == 500  & SAMPLE_NO == 31),
                                !(PROJECT_CODE == 'ABC' & SITE_ID == 505  & SAMPLE_NO == 33),
                                !(PROJECT_CODE == 'ABC' & SITE_ID == 505  & SAMPLE_NO == 36),
                                !(PROJECT_CODE == 'WDH' & SITE_ID == 9098 & SAMPLE_NO == 13))
} else if (attribute == 'col_p'){
  land_disturb <- c('5', '6', '7', '8')
  extract <- extract %>% filter(VALUE <= 150,
                                UD == 0,
                                between(LD, 9, 11), # LD == 10 - not working???
                                !grepl(paste(land_disturb, collapse='|'), DISTURB_TYPE),
                                !is.na(LATITUDE))
  extract <- extract %>% dplyr::select(LATITUDE, LONGITUDE, ID, VALUE)
  extract <- extract %>% dplyr::rename(X0.10.cm = VALUE)
}
# write.csv(extract, paste0(extd, '/lab_data.csv'), row.names=F)

#################################################################################################################
# load and read in lab data
lab.data <- extract
labdata1 <- reshape::rename(lab.data, c(UD = 'SUD', LD = 'SLD'))

##Remove sites with only one analysed depth
count1    <- count(labdata1, vars='ID')
morethan1 <- subset(count1, freq > 1)
labdata   <- merge(labdata1, morethan1, by=c('ID'))

## A inserts
Inserts       <- read.csv(fullfilename_A_inserts) #Import Anew_depths.csv data generated in SQLDev for each duplex soil 
Allvalue      <- merge(Inserts, labdata, by=c('PROJECT_CODE', 'SITE_ID')) #Add lab data to new_depths (one to one)
Allvalue$diff <- Allvalue$UD - Allvalue$SUD #Add a 'diff' field to Allvalue df which is the 'difference between sample depth and A/B horizon depth'
Value         <- subset(Allvalue, diff >= 0 & UD >= 10, select=c(ID, UD, LD, VALUE, diff)) #Limit records to lab samples in A horizon using 'diff' field and only where A horizon >= 10cm depth
Nearsamp      <- aggregate(x=Value$diff, by=list(ID=Value$ID), min) #Limit records to sample nearest to A/B change
Nearsamp1     <- reshape::rename(Nearsamp, c(x='diff')) #Rename column 'x' to 'diff' in 'Nearsamp' df to conform to 'Value' df
Aresult       <- join(Value, Nearsamp1, by=c('ID', 'diff'), type='right', match='first') #Combine df with nearest sample 'nearsamp' with df with lab data 'Value'
Aresult       <- subset(Aresult, select=c(ID, UD, LD, VALUE)) #Remove 'diff' in prep for Spline program

## B inserts
Inserts       <- read.csv(fullfilename_B_inserts) #Import Anew_depths.csv data generated in SQLDev for each duplex soil 
Allvalue      <- merge(Inserts, labdata, by=c('PROJECT_CODE', 'SITE_ID')) #Add lab data to new_depths (one to one)
Allvalue$diff <- Allvalue$LD - Allvalue$SLD #Add a 'diff' field to Allvalue df which is the 'difference between sample depth and A/B horizon depth'
Value         <- subset(Allvalue, diff <= 0 & UD >= 10, select=c(ID, UD, LD, VALUE, diff)) #Limit records to lab samples in B horizon using 'diff' field and only where B horizon Upper Depth >= 10cm
Nearsamp      <- aggregate(x=Value$diff, by=list(ID=Value$ID), max) #Limit records to sample nearest to A/B change
Nearsamp1     <- reshape::rename(Nearsamp, c(x='diff')) #Rename column 'x' to 'diff' in 'Nearsamp' df to conform to 'Value' df
Bresult       <- join(Value, Nearsamp1, by=c('ID', 'diff'), type='right', match='first') #Combine df with nearest sample 'nearsamp' with df with lab data 'Value'
Bresult       <- subset(Bresult, select=c(ID, UD, LD, VALUE)) #Remove 'diff' in prep for Spline program

## Real sample depths
Realdata    <- na.omit(subset(labdata, select=c(ID, SUD, SLD, VALUE))) ##select only fields required for spline and omit nulls 
DupsRemoved <- unique(Realdata[,1:2]) ## delete duplicates or overlaps
Realdata    <- join(Realdata, DupsRemoved, by=c('ID', 'SUD'), type='right', match='first')
Realdata    <- reshape::rename(Realdata, c(SUD='UD', SLD='LD')) #rename depth columns to conform with Aresult and Bresult df

##merge Real samples with A & B inserts
Result1 <- rbind(Realdata, Aresult, Bresult) #Append data from each result above into the one df
Result2 <- aggregate(x=Result1$VALUE, by=list(ID=Result1$ID, UD=Result1$UD, LD=Result1$LD), mean) #For identical sample depths with more than one measured value, take the average of these values
Result2 <- reshape::rename(Result2, c(x='VALUE'))
Result  <- Result2[order(Result2$ID, Result2$UD, Result2$LD),] #Order df records on ID, upper depth
# write.csv(Result, file=paste0(resd, attribute, '_with_inserts.csv'), row.names=F) #Export result to project directory

#################################################################################################################
# Spline fitting for horizon data (Matlab Code converted to R by Brendan Malone)
data <- Result
ea.fit.cs <- ea_spline(obj = data, var.name = "VALUE", d = t(c(0, 5, 15, 30, 60, 100, 200)), lam = 0.1, vlow = 0, show.progress = TRUE)
names(ea.fit.cs$harmonised) <- c('ID', 'X0.5.cm', 'X5.15.cm', 'X15.30.cm', 'X30.60.cm', 'X60.100.cm', 'X100.200.cm', 'soil.depth')
cords <- lab.data
cords <- subset(cords, select = c(ID, LATITUDE, LONGITUDE))
cords <- unique(cords)
fff <- merge(ea.fit.cs$harmonised, cords, by = c("ID"), all.y=T)
fff <- fff[c(9,10,1,2,3,4,5,6,7,8)]
fff <- reshape::rename(fff, c(LATITUDE = "Y", LONGITUDE = "X"))

data <- fff
# clean dataset
for (i in 4:ncol(data)){
  if (any(attribute == c('clay', 'cs', 'fs', 'silt', 'sand'))){
    data[i] <- replace(data[,i], !between(data[,i], 0, 100), NA)
  } else if (any(attribute == c('phw', 'phcl'))){
    data[i] <- replace(data[,i], !between(data[,i], 0, 14), NA)
  } else if (attribute == 'cat_ca'){
    data[i] <- replace(data[,i], !between(data[,i], 0, 200), NA)
  } else if (attribute == 'cat_cec'){
    data[i] <- replace(data[,i], !between(data[,i], 0, 300), NA)
  } else if (attribute == 'cat_k'){
    data[i] <- replace(data[,i], !between(data[,i], 0, 15), NA)
  } else if (attribute == 'cat_na'){
    data[i] <- replace(data[,i], !between(data[,i], 0, 150), NA)
  } else if (attribute == 'cat_mg'){
    data[i] <- replace(data[,i], !between(data[,i], 0, 100), NA)
  } else if (attribute == 'tn'){
    data[i] <- replace(data[,i], !between(data[,i], 0, 2), NA)
  } else if (attribute == 'tp'){
    data[i] <- replace(data[,i], !between(data[,i], 0, 2), NA)
  } else if (attribute == 'ec'){
    # data[i] <- replace(data[,i], !between(data[,i], 0, 100), NA)
  } else if (attribute == 'clay_act'){
    # data[i] <- replace(data[,i], !between(data[,i], 0, 100), NA)
  } else if (attribute == 'esp'){
    # data[i] <- replace(data[,i], !between(data[,i], 0, 100), NA)
  } else if (attribute == 'sar'){
    # data[i] <- replace(data[,i], !between(data[,i], 0, 100), NA)
  } else if (attribute == 'col_p'){
    data[i] <- replace(data[,i], !between(data[,i], 0, 150), NA)
  } else if (attribute == 'wb_oc'){
    # data[i] <- replace(data[,i], !between(data[,i], 0, 100), NA)
  }
}

# if (attribute == 'col_p'){
#   qlump <- vect('/scratch/rsc3/leos/QSRS/QLD_LANDUSE_June_2019/Current_landuse.shp')
#   codes <- read.csv('/scratch/rsc3/leos/QSRS/QLD_LANDUSE_June_2019/ALUMClassV8.csv')
#   qlump <- merge(qlump, codes, by.x='ALUM_Code', by.y='Tertiary')
#   qlump <- subset(qlump, qlump$FertYN == 2)
#   qlump_buf <- terra::buffer(qlump, width=75)
#   col_p <- st_as_sf(data, coords=c('LONGITUDE','LATITUDE'), crs=crs('EPSG:4326')) %>% st_transform(crs=crs(qlump)) %>% vect
#   SitesExcluded <- col_p[qlump_buf,]
#   ExcludedSites <- as.data.frame(SitesExcluded) 
#   ExcludedSites <- subset(ExcludedSites, X0.10.cm > 7)
#   FinalSites <- col_p[!(col_p$ID %in% ExcludedSites$ID),]
#   FinalSites <- terra::project(FinalSites, 'EPSG:4326')
#   final_pts <- as.data.frame(FinalSites)
#   coords <- crds(FinalSites)
#   data <- data.frame('LATITUDE' = coords[,2], 'LONGITUDE' = coords[,1], 'ID' = final_pts$ID, 'X0.10.cm' = final_pts$X0.10.cm)
# }

write.csv(data, paste0(extd, '/harmonized.csv'), row.names=F)
