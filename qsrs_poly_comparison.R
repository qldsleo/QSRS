# load in qsrs polygon
qsrs_poly <- vect(paste0(code,'BEST_AVAIL_SOIL_POLY_041022.gdb'))
qsrs_poly <- qsrs_poly[1:30,]

# load in qsrs rasters
qsrs_files <- list.files(paste0(resd, 'maps/', ifelse(run_PCA,'PCA','all')), full.names=T, pattern='mean')
qsrs_rast <- rast(qsrs_files)
names(qsrs_rast) <- gsub('.tif', '', list.files(paste0(resd, 'maps/', ifelse(run_PCA,'PCA','all')), pattern='mean'))
qsrs_rast <- terra::crop(qsrs_rast, qsrs_poly, mask=T)
plot(qsrs_rast)

qsrs_rast_extract <- terra::extract(qsrs_rast, qsrs_poly, fun=mean, bind=T, na.rm=T)
extract_df <- qsrs_rast_extract %>% as.data.frame %>% mutate(
  CS1_UPPERDEPTH = CS1_UPPERDEPTH * 100,
  CS1_LOWERDEPTH = CS1_LOWERDEPTH * 100,
  CS2_UPPERDEPTH = CS2_UPPERDEPTH * 100,
  CS2_LOWERDEPTH = CS2_LOWERDEPTH * 100,
  CS3_UPPERDEPTH = CS3_UPPERDEPTH * 100,
  CS3_LOWERDEPTH = CS3_LOWERDEPTH * 100,
  CS4_UPPERDEPTH = CS4_UPPERDEPTH * 100,
  CS4_LOWERDEPTH = CS4_LOWERDEPTH * 100,
  CS5_UPPERDEPTH = CS5_UPPERDEPTH * 100,
  CS5_LOWERDEPTH = CS5_LOWERDEPTH * 100
)

test <- extract_df %>% mutate(
  CS1_CLAY_QSRS = case_when(CS1_UPPERDEPTH == 0 & CS1_LOWERDEPTH <= 9  ~ clay_mean_0.5cm_30m_v2,
                            CS1_UPPERDEPTH == 0 & between(CS1_LOWERDEPTH,6,22) ~ mean(c(clay_mean_0.5cm_30m_v2,clay_mean_5.15cm_30m_v2))),
  CS3_CLAY_QSRS = case_when()
)

ggplot(test, aes(CS1_CLAY, CS1_CLAY_QSRS))+geom_point()+geom_abline(colour='red')+
  scale_y_continuous(limits=c(0,60))+scale_x_continuous(limits=c(0,60))



extract_df$CS3_UPPERDEPTH
extract_df$CS3_LOWERDEPTH
extract_df$CS3_CLAY


test$CS1_CLAY_QSRS


qsrs_poly_data <- qsrs_poly %>% filter()
qsrs_poly$CS1_LOWERDEPTH