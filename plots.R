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