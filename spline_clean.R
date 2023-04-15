extract <- read.csv(fullfilename_td)
extract <- extract %>% select(ID, UD, LD, paste0(attribute), paste0(attribute, '_method'), OBJECTID, PROJECT_CODE, SITE_ID, OBS_NO, HORIZON_NO, SAMPLE_NO,
                              YEAR, LATITUDE, LONGITUDE, ZONE, EASTING, NORTHING, DISTURB_TYPE, NATURE)
names(extract)[which(names(extract) == paste0(attribute))] <- "VALUE"
extract <- extract %>% filter(!is.na(VALUE))
all_extract <- distinct(extract, ID, UD, LD, .keep_all=T)

#################################################################################################################
# load and read in lab data
lab.data <- all_extract
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
write.csv(Result, file=paste0(resd, attribute, '_with_inserts.csv'), row.names=F) #Export result to project directory

#################################################################################################################

# Spline fitting for horizon data (Matlab Code converted to R by Brendan Malone)
ea_spline <- function(obj, var.name, lam = 0.1, d = t(c(0,5,15,30,60,100,200)), vlow = 0, vhigh = 1000, show.progress=TRUE){
  if (is(obj,"SoilProfileCollection") == TRUE){
    depthcols = obj@depthcols
    idcol = obj@idcol
    obj@horizons = obj@horizons[,c(idcol, depthcols, var.name)]
    d.dat<- as.data.frame(obj@horizons)}
  
  if (is(obj,"data.frame") == TRUE){
    d.dat<- as.data.frame(obj[,c(1:3,which(colnames(obj)==var.name))])}
  
  if (is(obj,"data.frame") == FALSE & is(obj,"SoilProfileCollection") == FALSE){
    stop("ERROR:: Data must be either class data.frame or SoilProfileCollection")}
  
  mxd<- max(d)
  s2<- (0.05*mean(d.dat[,4]))^2  # overall variance of soil attribute
  sp_dat<-split(d.dat,d.dat[,1]) 
  nl<- length(lam)  # Length of the lam matrix
  
  # matrix of the continous splines for each data point
  m_fyfit<- matrix(NA,ncol=length(c(1:mxd)),nrow=length(sp_dat))
  
  # Matrix in which the averaged values of the spline are fitted. The depths are specified in the (d) object
  yave<- matrix(NA,ncol=length(d),nrow=length(sp_dat))
  
  # Matrix in which the sum of square errors of each lamda iteration for the working profile are stored
  sse<- matrix(NA,ncol=length(lam),nrow=1)
  
  # Matrix in which the sum of square errors for eac h lambda iteration for each profile are stored
  sset<- matrix(NA,ncol=length(lam),nrow=length(sp_dat))
  
  #Profile ids
  mat_id<- d.dat[0,]
  
  #combined data frame for observations and spline predictions
  dave<- d.dat[1,]
  dave$predicted<- 0
  dave$FID<- 0
  dave<- dave[0,]
  
  
  
  ## Fit splines profile by profile:            
  message("Fitting mass preserving splines per profile...")
  if (show.progress) pb <- txtProgressBar(min=0, max=length(sp_dat), style=3)
  cnt<- 1
  for(st in 1:length(sp_dat)) {
    subs<-sp_dat[[st]]  # subset the profile required
    subs<-as.matrix(subs)
    mat_id[st,1]<- subs[1,1]
    
    
    # manipulate the profile data to the required form
    ir<- c(1:nrow(subs))
    ir<-as.matrix(t(ir))
    u<- as.numeric(subs[ir,2])
    u<-as.matrix(t(u))   # upper 
    v<- as.numeric(subs[ir,3])
    v<-as.matrix(t(v))   # lower
    y<- as.numeric(subs[ir,4])
    y<-as.matrix(t(y))   # concentration 
    n<- length(y);       # number of observations in the profile
    
    
    ############################################################################################################################################################
    ### routine for handling profiles with one observation
    if (n == 1)
    {dave[cnt:(cnt-1+nrow(subs)),1:4]<- subs
    dave[cnt:(cnt-1+nrow(subs)),5]<- y
    dave[cnt:(cnt-1+nrow(subs)),6]<- st
    xfit<- as.matrix(t(c(1:mxd))) # spline will be interpolated onto these depths (1cm res)
    nj<- max(v)
    if (nj > mxd)
    {nj<- mxd}
    yfit<- xfit
    yfit[,1:nj]<- y   # values extrapolated onto yfit
    if (nj < mxd)
    {yfit[,(nj+1):mxd]=-9999}
    m_fyfit[st,]<- yfit
    
    
    # Averages of the spline at specified depths
    nd<- length(d)-1  # number of depth intervals
    dl<-d+1     #  increase d by 1
    for (cj in 1:nd) { 
      xd1<- dl[cj]
      xd2<- dl[cj+1]-1
      if (nj>=xd1 & nj<=xd2)
      {xd2<- nj-1
      yave[st,cj]<- mean(yfit[,xd1:xd2])}
      else
      {yave[st,cj]<- mean(yfit[,xd1:xd2])}   # average of the spline at the specified depth intervals
      yave[st,cj+1]<- max(v)} #maximum soil depth
    cnt<- cnt+nrow(subs)
    ##################################
    }
    # End of single observation profile routine
    ###############################################################################################################################################################
    
    ## Start of routine for fitting spline to profiles with multiple observations         
    
    else  {    
      ###############################################################################################################################################################
      dave[cnt:(cnt-1+nrow(subs)),1:4]<- subs
      dave[cnt:(cnt-1+nrow(subs)),6]<- st
      ## ESTIMATION OF SPLINE PARAMETERS
      np1 <- n+1  # number of interval boundaries
      nm1 <- n-1
      delta <- v-u  # depths of each layer
      del <- c(u[2:n],u[n])-v   # del is (u1-v0,u2-v1, ...)
      
      ## create the (n-1)x(n-1) matrix r; first create r with 1's on the diagonal and upper diagonal, and 0's elsewhere
      r <- matrix(0,ncol=nm1,nrow=nm1)
      for(dig in 1:nm1){
        r[dig,dig]<-1
      }
      for(udig in 1:nm1-1){
        r[udig,udig+1]<-1
      }
      
      ## then create a diagonal matrix d2 of differences to premultiply the current r
      d2 <- matrix(0, ncol=nm1, nrow=nm1)
      diag(d2) <- delta[2:n]  # delta = depth of each layer
      
      ## then premultiply and add the transpose; this gives half of r
      r <- d2 %*% r
      r <- r + t(r)
      
      ## then create a new diagonal matrix for differences to add to the diagonal
      d1 <- matrix(0, ncol=nm1, nrow=nm1)
      diag(d1) <- delta[1:nm1]  # delta = depth of each layer
      
      d3 <- matrix(0, ncol=nm1, nrow=nm1)
      diag(d3) <- del[1:nm1]  # del =  differences
      
      r <- r+2*d1 + 6*d3
      
      ## create the (n-1)xn matrix q
      q <- matrix(0,ncol=n,nrow=n)
      for (dig in 1:n){
        q[dig,dig]<- -1 
      }
      for (udig in 1:n-1){
        q[udig,udig+1]<-1 
      }
      q <- q[1:nm1,1:n]
      dim.mat <- matrix(q[],ncol=length(1:n),nrow=length(1:nm1))
      
      ## inverse of r
      rinv <- try(solve(r), TRUE)
      
      ## Note: in same cases this will fail due to singular matrix problems, hence you need to check if the object is meaningfull:
      if(is.matrix(rinv)){
        ## identity matrix i
        ind <- diag(n)
        
        ## create the matrix coefficent z
        pr.mat <- matrix(0,ncol=length(1:nm1),nrow=length(1:n))
        pr.mat[] <- 6*n*lam
        fdub <- pr.mat*t(dim.mat)%*%rinv
        z <- fdub%*%dim.mat+ind
        
        ## solve for the fitted layer means
        sbar <- solve(z,t(y))
        dave[cnt:(cnt-1+nrow(subs)),5]<- sbar
        cnt<- cnt+nrow(subs)
        
        
        ## calculate the fitted value at the knots
        b <- 6*rinv%*%dim.mat%*% sbar
        b0 <- rbind(0,b) # add a row to top = 0
        b1 <- rbind(b,0) # add a row to bottom = 0
        gamma <- (b1-b0) / t(2*delta)
        alfa <- sbar-b0 * t(delta) / 2-gamma * t(delta)^2/3
        
        ## END ESTIMATION OF SPLINE PARAMETERS
        ###############################################################################################################################################################
        
        
        ## fit the spline 
        xfit<- as.matrix(t(c(1:mxd))) ## spline will be interpolated onto these depths (1cm res)
        nj<- max(v)
        if (nj > mxd)
        {nj<- mxd}
        yfit<- xfit
        for (k in 1:nj){
          xd<-xfit[k]
          if (xd< u[1])
          {p<- alfa[1]} else
          {for (its in 1:n) {
            if(its < n)
            {tf2=as.numeric(xd>v[its] & xd<u[its+1])} else {tf2<-0}
            if (xd>=u[its] & xd<=v[its])
            {p=alfa[its]+b0[its]*(xd-u[its])+gamma[its]*(xd-u[its])^2} else if (tf2)
            {phi=alfa[its+1]-b1[its]*(u[its+1]-v[its])
            p=phi+b1[its]*(xd-v[its])}
          }}
          yfit[k]=p }
        if (nj < mxd)
        {yfit[,(nj+1):mxd]=NA}
        
        yfit[which(yfit > vhigh)] <- vhigh
        yfit[which(yfit < vlow)]  <-vlow
        m_fyfit[st,]<- yfit
        
        ## Averages of the spline at specified depths
        nd<- length(d)-1  # number of depth intervals
        dl<-d+1     #  increase d by 1
        for (cj in 1:nd) { 
          xd1<- dl[cj]
          xd2<- dl[cj+1]-1
          if (nj>=xd1 & nj<=xd2)
          {xd2<- nj-1
          yave[st,cj]<- mean(yfit[,xd1:xd2])}
          else
          {yave[st,cj]<- mean(yfit[,xd1:xd2])}   # average of the spline at the specified depth intervals
          yave[st,cj+1]<- max(v)} #maximum soil depth 
        
        
        
        
        ## CALCULATION OF THE ERROR BETWEEN OBSERVED AND FITTED VALUES
        ## calculate Wahba's estimate of the residual variance sigma^2
        ssq <- sum((t(y)-sbar)^2)
        g <- solve(z)
        ei <- eigen(g)
        ei <- ei$values
        df <- n-sum(ei)
        sig2w <- ssq/df
        ## calculate the Carter and Eagleson estimate of residual variance
        dfc <- n-2*sum(ei)+sum(ei^2)
        sig2c <- ssq/dfc
        ## calculate the estimate of the true mean squared error
        tmse <- ssq/n-2*s2*df/n+s2
        sset[st] <- tmse
        
      }
    }
    
    if (show.progress) { setTxtProgressBar(pb, st)  }
  }
  if (show.progress) { 
    close(pb)
    
  }
  
  
  ## asthetics for output 
  ## yave
  yave<- as.data.frame(yave)
  jmat<- matrix(NA,ncol=1,nrow=length(d))
  for (i in 1:length(d)-1) {
    a1<-paste(d[i],d[i+1],sep="-")
    a1<-paste(a1,"cm",sep=" ")
    jmat[i]<- a1}
  jmat[length(d)]<- "soil depth"
  for (jj in 1:length(jmat)){
    names(yave)[jj]<- jmat[jj] 
  }
  yave<- cbind(mat_id[,1],yave)
  names(yave)[1]<- "id"
  
  retval <- list(harmonised=yave, obs.preds=dave,var.1cm=t(m_fyfit),tmse=sset)
  
  return(retval)
  
  
}

# load data
data <- Result
ea.fit.cs <- ea_spline(obj = data, var.name = "VALUE", d = t(c(0, 5, 15, 30, 60, 100, 200)), lam = 0.1, vlow = 0, show.progress = TRUE)
names(ea.fit.cs$harmonised) <- c('ID', 'X0.5.cm', 'X5.15.cm', 'X15.30.cm', 'X30.60.cm', 'X60.100.cm', 'X100.200.cm', 'soil.depth')
cords <- lab.data
cords <- subset(cords, select = c(ID, LATITUDE, LONGITUDE))
cords <- unique(cords)
fff <- merge(ea.fit.cs$harmonised, cords, by = c("ID"), all.y=T)
fff <- fff[c(9,10,1,2,3,4,5,6,7,8)]
fff <- reshape::rename(fff, c(LATITUDE = "Y", LONGITUDE = "X"))
write.csv(fff, paste0(save_location, 'harmonized.csv'), row.names=F)


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

if (attribute == 'col_p'){
  qlump <- vect('/scratch/rsc3/leos/QSRS/QLD_LANDUSE_June_2019/Current_landuse.shp')
  codes <- read.csv('/scratch/rsc3/leos/QSRS/QLD_LANDUSE_June_2019/ALUMClassV8.csv')
  qlump <- merge(qlump, codes, by.x='ALUM_Code', by.y='Tertiary')
  qlump <- subset(qlump, qlump$FertYN == 2)
  qlump_buf <- terra::buffer(qlump, width=75)
  col_p <- st_as_sf(data, coords=c('LONGITUDE','LATITUDE'), crs=crs('EPSG:4326')) %>% st_transform(crs=crs(qlump)) %>% vect
  SitesExcluded <- col_p[qlump_buf,]
  ExcludedSites <- as.data.frame(SitesExcluded) 
  ExcludedSites <- subset(ExcludedSites, X0.10.cm > 7)
  FinalSites <- col_p[!(col_p$ID %in% ExcludedSites$ID),]
  FinalSites <- terra::project(FinalSites, 'EPSG:4326')
  final_pts <- as.data.frame(FinalSites)
  coords <- crds(FinalSites)
  data <- data.frame('LATITUDE' = coords[,2], 'LONGITUDE' = coords[,1], 'ID' = final_pts$ID, 'X0.10.cm' = final_pts$X0.10.cm)
}