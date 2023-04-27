range01 <- function(x){
  (x-min(x))/(max(x)-min(x))
}

d = 1

pred_depth <- depths[d]
model <- readRDS(paste0(resd, method_dsm, '/models/', pred_depth, '/mapMod.rds'))
model

DSM_data <- read.csv(paste0(extd, '/Site_Covariate_intersect.csv'))
DSM_data <- DSM_data %>% drop_na(names(DSM_data)[10:ncol(DSM_data)])
if(!is.null(covs_subset_cat)){
  for(fname_cov in covs_subset_cat){
    covname_Tmp <- gsub('.tif', '', fname_cov)
    if (covname_Tmp %in% names(DSM_data)){
      DSM_data[,covname_Tmp] <- as.factor(DSM_data[,covname_Tmp])
      print(paste0('Converted - ', covname_Tmp, ' - to factor'))
    }
  }
}

shap <- fastshap::explain(model, X = DSM_data[10:ncol(DSM_data)], nsim = 1, pred_wrapper = predict)
shap

autoplot(shap)
autoplot(shap, type = "contribution", row_num = 1)
autoplot(shap, type = "dependence", feature = "PM_radmap_v4_2019_filtered_dose_GAPFilled_qld", X = DSM_data, smooth = T, smooth_color = 'blue')

shap2 <- shap %>% as.data.frame %>% mutate(ID = 1:nrow(shap))
shap_pivot <- shap2 %>% pivot_longer(!ID, names_to = 'covariate', values_to = 'shapley_values')

DSM_data <- read.csv(paste0(extd, '/Site_Covariate_intersect.csv'))
DSM_data <- DSM_data %>% drop_na(names(DSM_data)[10:ncol(DSM_data)])
shap_norm <- DSM_data[10:ncol(DSM_data)] %>% mutate_all(range01)
if(!is.null(covs_subset_cat)){
  for(fname_cov in covs_subset_cat){
    covname_Tmp <- gsub('.tif', '', fname_cov)
    if (covname_Tmp %in% names(shap_norm)){
      shap_norm[,covname_Tmp] <- NA
      print(paste0('Converted - ', covname_Tmp, ' - to factor'))
    }
  }
}
shap_norm <- shap_norm %>% mutate(ID = 1:nrow(shap))
shap_pivot_norm <- shap_norm %>% pivot_longer(!ID, names_to = 'covariate', values_to = 'shapley_values_norm')

shap_pivot$shapley_values_norm <- shap_pivot_norm$shapley_values_norm
ggplot(shap_pivot)+geom_point(aes(y=covariate, x=shapley_values, colour=shapley_values_norm))+
  theme_bw()+
  labs(x='Shapley Value', y=NULL)+
  scale_colour_viridis_c(option = "plasma", direction = -1)+
  geom_vline(xintercept = 0)
