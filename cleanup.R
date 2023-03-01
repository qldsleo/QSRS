# clean up tile folders from mapping script

# for (d in 1:length(depths)){
for (d in 1:6){
    
  # which depth
  pred_depth <- depths[d]
  
  unlink(paste0(resd, 'maps/', ifelse(run_PCA,'PCA/','all/'), pred_depth), recursive=T)

}
