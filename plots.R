plot_data <- do.call(rbind, lapply(1:length(depths), function(i){
  read.csv(paste0(save_dir, '/', depths[i], '/cDat_GOOF.csv'))
}))
plot_data$Depth = factor(plot_data$Depth, levels=c(depths))
ggplot(plot_data, aes(Obs, Pred))+
  geom_point(size=0.5)+geom_abline(col='red', size=0.75)+geom_smooth(method='lm',formula='y ~ x', col='blue', size=0.75)+
  facet_wrap(~ Depth, scales='free')+
  theme_bw()
ggsave(filename=paste0(save_dir, '/', 'Cal_alldepths.jpeg'))

plot_data <- do.call(rbind, lapply(1:length(depths), function(i){
  read.csv(paste0(save_dir, '/', depths[i], '/vDat_GOOF.csv'))
}))
plot_data$Depth = factor(plot_data$Depth, levels=c(depths))
ggplot(plot_data, aes(Obs, Pred))+
  geom_point(size=0.5)+geom_abline(col='red', size=0.75)+geom_smooth(method='lm',formula='y ~ x', col='blue', size=0.75)+
  facet_wrap(~ Depth, scales='free')+
  theme_bw()
ggsave(filename=paste0(save_dir, '/','Val_alldepths.jpeg'))

# PICP plots
GOOFDat_long <- GOOFDat[1:length(depths),] %>% 
  pivot_longer(-!c(PICP_5, PICP_10, PICP_20, PICP_40, PICP_60, PICP_80, PICP_90, PICP_95, PICP_97.5, PICP_99), names_to='PICP', values_to='picpVec') %>% as.data.frame
GOOFDat_long$piclplot <- piclplot
GOOFDat_long$Depths = factor(GOOFDat_long$Depths, levels=c(depths))

ggplot(GOOFDat_long, aes(piclplot, picpVec))+
  geom_point(size=0.5)+geom_abline(col='red', size=0.75)+
  facet_wrap(~ Depths, scales='free')+
  theme_bw()+labs(y='PICP', x='Confidence level')
ggsave(filename=paste0(save_dir, '/', 'PICP_alldepths.jpeg'))