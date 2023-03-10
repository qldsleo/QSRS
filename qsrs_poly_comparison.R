# load in qsrs polygon
qsrs_poly <- vect(paste0(covd,'/QSRS_poly.shp'))
qsrs_poly <- qsrs_poly[1:30,]

# load in qsrs rasters
qsrs_files <- list.files(paste0(resd, method_dsm, '/maps/', ifelse(run_PCA,'PCA','all')), full.names=T, pattern='mean')
qsrs_rast <- rast(qsrs_files)
names(qsrs_rast) <- gsub('.tif', '', list.files(paste0(resd, method_dsm, '/maps/', ifelse(run_PCA,'PCA','all')), pattern='mean'))
qsrs_rast <- terra::crop(qsrs_rast, qsrs_poly, mask=T)
plot(qsrs_rast)

qsrs_rast_extract <- terra::extract(qsrs_rast, qsrs_poly, fun=mean, bind=T, na.rm=T)
extract_df <- qsrs_rast_extract %>% as.data.frame %>% mutate(
  CS1_UPPERD = CS1_UPPERD * 100,
  CS1_LOWERD = CS1_LOWERD * 100,
  CS2_UPPERD = CS2_UPPERD * 100,
  CS2_LOWERD = CS2_LOWERD * 100,
  CS3_UPPERD = CS3_UPPERD * 100,
  CS3_LOWERD = CS3_LOWERD * 100,
  CS4_UPPERD = CS4_UPPERD * 100,
  CS4_LOWERD = CS4_LOWERD * 100,
  CS5_UPPERD = CS5_UPPERD * 100,
  CS5_LOWERD = CS5_LOWERD * 100
)

test <- extract_df %>% mutate(
  CS1_CLAY_QSRS = case_when(CS1_UPPERD == 0 & CS1_LOWERD <= 9  ~ clay_mean_0.5cm_30m_v3)
)

ggplot(test, aes(CS1_CLAY, CS1_CLAY_QSRS))+geom_point()+geom_abline(colour='red')+
  scale_y_continuous(limits=c(0,60))+scale_x_continuous(limits=c(0,60))+
  labs(x='QSRS clay polygon', y='QSRS clay raster')
