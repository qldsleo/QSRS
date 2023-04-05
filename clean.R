# clean dataset
for (i in 4:ncol(data)){
  if (any(attribute == c('clay', 'cs', 'fs', 'silt'))){
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