ccc <- function(observed, predicted){
  mx=mean(observed)
  my=mean(predicted)
  s2x=var(observed)
  s2y=var(predicted)
  sxy=mean((observed-mx)*(predicted-my))
  ccc=2*sxy/(s2x+s2y+(mx-my)^2 )
  return(ccc)
} # concordance function (M. Webb 2015)
summaryStatsFn <- function(observed, predicted){
  
  # Note, R2 in defaultSummary is cor(observed , predicted) ^ 2, not 1-sum((observed - predicted) ^ 2) / sum((observed - mean(observed)) ^ 2) 
  # this function uses the cor(observed , predicted) ^ 2 version of R2
  
  pvoVec <- rep(NA , 9)
  names(pvoVec) <- c('RMSE' , 'R2' , 'MAE' , 'MSE' , 'CC' , 'Bias' , 'RMedSE' , 'MedAE' , 'MedBias')
  
  pvoVec[1] <- sqrt(mean((observed - predicted) ^ 2)) #RMSE
  pvoVec[2] <- cor(observed , predicted) ^ 2 # R2 as sqr of cor (as given in defaultSummary)
  pvoVec[3] <- mean(abs(observed - predicted)) #MAE
  
  pvoVec[4] <- mean((observed - predicted)^2) #MSE
  pvoVec[5] <- ccc(observed, predicted) # Concordance
  pvoVec[6] <- sum(predicted - observed)/length(observed) # Bias
  pvoVec[7] <- sqrt(median((predicted - observed)^2)) # RMedSE
  pvoVec[8] <- median(abs(predicted - observed)) # MedAE
  pvoVec[9] <- median(predicted - observed) # MedBias
  
  return(pvoVec)
} # summary stats function
getTrainAndTestSamples_BalancedFactors <- function(data, percentTrain = propn_training, inSequence = F, seed = 0){
  
  set.seed(seed) # Set Seed so that same sample can be reproduced in future also
  sample <- NULL
  train <- NULL
  test <- NULL
  
  if (d == 7){
    data <- data %>% filter(Depth == 2.5)
  }
  
  listOfFactorsAndTheirLevels <- lapply(Filter(is.factor, data), summary)
  factorContainingOneElement <- sapply(listOfFactorsAndTheirLevels, function(x){any(x==1)})
  
  if (any(factorContainingOneElement)){
    factors <- as.data.frame(factorContainingOneElement)
    factors_to_drop <- factors %>% filter(factorContainingOneElement == T) %>% row.names()
    data[,c(factors_to_drop)] <- NULL
    print('Dropping covariates - ')
    print(factors_to_drop)
  }
  
  if (length(Filter(is.factor, data)) == 0){
    sample <- sample(nrow(DSM_data_depth), propn_training * nrow(DSM_data_depth)) # Vector with training data row numbers
    
  }else{
    
    # Repeat loop until all Factor variables have same levels on both train and test set
    repeat{
      # Now Selecting 'percentTrain' of data as sample from total 'n' rows of the data  
      if (inSequence)
        sample <- 1:floor(percentTrain * nrow(data))
      else
        sample <- sample.int(n = nrow(data), size = floor(percentTrain * nrow(data)), replace = F)    
      
      train <- data[sample, ]   # create train set
      test  <- data[-sample, ]  # create test set
      
      train_factor_only <- Filter(is.factor, train) # df containing only 'train' Factors as columns
      test_factor_only <- Filter(is.factor, test)   # df containing only 'test' Factors as columns
      
      haveFactorsWithExistingLevels <- NULL
      for (i in 1:ncol(train_factor_only)){ # for each column (i.e. factor variable)
        names_train <- names(summary(train_factor_only[,i]))  # get names of all existing levels in this factor column for 'train'
        names_test <- names(summary(test_factor_only[,i]))    # get names of all existing levels in this factor column for 'test'
        
        symmetric_diff <- union(setdiff(names_train,names_test), setdiff(names_test,names_train)) # get the symmetric difference between the two factor columns (from 'train' and 'test')
        if (length(symmetric_diff) == 0)  # if no elements in the symmetric difference set then it means that both have the same levels present at least once
          haveFactorsWithExistingLevels <- c(haveFactorsWithExistingLevels, TRUE) # append logic TRUE
        else # if some elements in the symmetric difference set then it means that one of the two sets (train, test) has levels the other doesn't and it will eventually flag up when using function predict()
          haveFactorsWithExistingLevels <- c(haveFactorsWithExistingLevels, FALSE) # append logic FALSE
      }
      
      if(all(haveFactorsWithExistingLevels))
        break # break out of the repeat loop because we found a split that has factor levels existing in both 'train' and 'test' sets, for all Factor variables
    }
  }
  
  # return (list(sample, factors_to_drop))
  if(!exists('factors_to_drop')){
    return(list(sample, NA))
  }else{
    return(list(sample, factors_to_drop))
  }
}