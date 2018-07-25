# some functions shared by the analyses of the two types of data

# localization data can be downloaded from OSF:
nocursorURLs <- c('exposure'='https://osf.io/9s6au/?action=download', 'classic'='https://osf.io/8hm7f/?action=download')

# no-cursor data can be downloaded from OSF:
localizationURLs <- c('exposure'='https://osf.io/bk4s6/?action=download', 'classic'='https://osf.io/upw49/?action=download', 'online'='https://osf.io/wjcgk/download')

# information on exposure participants
informationURLs <- c('blinkdetect'='https://osf.io/yrpkm/?action=download', 'demographics'='https://osf.io/7sn4b/?action=download')

# we'll control the color in figures centrally from here:
colorset <- list()

colorset[['expActS']] <- '#005de4ff' # blue
colorset[['expActT']] <- '#005de42f'
colorset[['expPasS']] <- '#0fd2e2ff' # lighter blue
colorset[['expPasT']] <- '#0fd2e22f'
colorset[['claActS']] <- '#e51636ff' # "York red"
colorset[['claActT']] <- '#e516362f'
colorset[['claPasS']] <- '#ff8200ff' # orange
colorset[['claPasT']] <- '#ff82002f'

# colorset[['extra1S']] <- '#c400c4ff' # purple
# colorset[['extra1T']] <- '#c400c42f'

colorset[['onlActS']] <- '#b400e4ff' # purple
colorset[['onlActT']] <- '#b400e42f'
# colorset[['onlPasS']] <- '#8266f4ff' # violet
# colorset[['onlPasT']] <- '#8266ff2f'
colorset[['onlPasS']] <- '#ff6ec7ff' # pink
colorset[['onlPasT']] <- '#ff6ec72f'

library('ez')

installRequire.Packages <- function(packages) {
  
  installed.list <- rownames(installed.packages())
  missingpackages <- c()
  
  for (pkg in packages) {
    
    if (!pkg %in% installed.list) {
      # we can't force people to install packages, instead, we'll warn them
      #install.packages(pkg,dep=TRUE)
      missingpackages <- c(missingpackages, pkg)
    }
    
    require(pkg, character.only=TRUE)
    
  }
  
  if (length(missingpackages) > 0) {
    cat('\nWARNING, some pacakages are missing:\n')
    cat(missingpackages)
    cat('\nlme4 and lmerTest are required to reproduce the Satterthwaite LME analyses\n')
    cat('\ncar and nlme are required to reproduce the analyses using Chi-squared apprcimation\n')
    cat('\nforeach and doParallel speed up a few functions\n')
    cat('\nsvglite is used to produce SVG files with figures\n')
  }
  
}

load.DownloadDataframe <- function(url,filename) {
  
  if (file.exists(filename)) {
    
    df <- read.csv(filename, stringsAsFactors=FALSE)
    
  } else {
    
    df <- read.csv(url(url),stringsAsFactors=FALSE)
    
    write.csv(df,filename,row.names=FALSE,quote=FALSE)
    
  }
  
  return(df)
  
}

t.interval = function(data, variance = var(data, na.rm=TRUE), conf.level = 0.95) {
  
  data <- data[!is.na(data)]
  
  z = qt((1 - conf.level)/2, df = length(data) - 1, lower.tail = FALSE)
  
  xbar = mean(data)
  sdx = sqrt(variance/length(data))
  
  return(c(xbar - z * sdx, xbar + z * sdx))
  
}


BS.interval = function(data, conf.level = .95, bootstraps = 1000) {
  
  resampled <- matrix(data=sample(data,size=length(data)*bootstraps,replace=TRUE),ncol=bootstraps)
  
  BSmeans <- apply(resampled,MARGIN=2,FUN=mean)
  
  oneside <- (1 - conf.level) / 2
  return(quantile(BSmeans,probs=c(oneside,1-oneside)))
  
}

bootstrapGaussianPeak <- function(data,bootstraps=1000,mu=47.5,sigma=30,scale=10,offset=4,CIs=c(.95)) {
  
  # parallel for-loop?
  # installRequire.Packages(c('foreach','doParallel'))
  
  installed.list <- rownames(installed.packages())
  packages <- c('foreach','doParallel')
  allInstalled <- TRUE
  for (pkg in packages) {
    if (!pkg %in% installed.list) {
      allInstalled <- FALSE
    } else {
      require(pkg, character.only=TRUE)
    }
  }
  
  # these are.... the hand angles / targets?
  x <- as.numeric(names(colMeans(data)))
  
  if (allInstalled) {
    cores=detectCores()
    
    # very friendly for the rest of the system
    usecores <- max(1, ceiling((cores[1] / 2) - 1))
    # somewhat unfriendly for the rest of the system
    usecores <- max(1, (cores[1] - 1))
    
    cl <- makeCluster(usecores)
    registerDoParallel(cl)
    
    mus <- foreach (iteration=1:bootstraps, .combine=rbind, .export=c('getGaussianFit','GaussianErrors','parGaussian')) %dopar% {
      
      y <- colMeans(data[sample(row.names(data),size=nrow(data),replace=TRUE),])
      
      as.numeric(getGaussianFit(x,y,mu=47.5,sigma=30,scale=10,offset=4)$par['mu'])
      
    }
    
    stopCluster(cl)
    
  } else {
    
    mus <- c()
    
    for (iteration in c(1:bootstraps)) {
      
      y <- colMeans(data[sample(row.names(data),size=nrow(data),replace=TRUE),])
      
      mus <- c(mus, as.numeric(getGaussianFit(x,y,mu=47.5,sigma=30,scale=10,offset=4)$par['mu']))
      
    }
    
  }
  
  probs <- c(.50)
  for (CI in CIs) {
    probs <- c(((1-CI)/2), probs, 1-((1-CI)/2))
  }
  # c(.025,.05,.10,.50,.90,.95,.975)
  return(quantile(mus,probs=probs))
  
}

# read this vignette:
# https://cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf
# 
# use this code from SO:
# 
# library(foreach)
# library(doParallel)
# 
# #setup parallel backend to use many processors
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
# 
# finalMatrix <- foreach(i=1:150000, .combine=cbind) %dopar% {
#   tempMatrix = functionThatDoesSomething() #calling a function
#   #do other things if you want
#   
#   tempMatrix #Equivalent to finalMatrix = cbind(finalMatrix, tempMatrix)
# }
# #stop cluster
# stopCluster(cl)


getGaussianFit <- function(x,y,mu=47.5,sigma=30,scale=10,offset=4) {
  
  data = data.frame(x,y)
  
  optim(par=c('mu'=mu,'sigma'=sigma, 'scale'=scale, 'offset'=offset), GaussianErrors, data=data)
  
}

GaussianErrors <- function(par, data) {
  
  errors <- sum(abs(data$y - (par['offset']+(par['scale']*parGaussian(par, data$x)))))^2
  
  return(errors)
  
}

parGaussian <- function(par,x) {
  
  y <- (1/(par['sigma']*sqrt(2*pi)))*exp(-0.5*(((x-par['mu'])/par['sigma']))^2)
  
  return(y)
  
}


# handling localization data -----
#   exp <- getPointLocalization('exposure', difference=TRUE, verbose=FALSE, selectPerformance=selectPerformance)

getPointLocalization <- function(group, difference=TRUE, points=c(15,25,35,45,55,65,75), movementtype='both', verbose=TRUE, LRpart='all', selectPerformance=selectPerformance) {
  
  df <- load.DownloadDataframe(url=localizationURLs[group],filename=sprintf('localization_%s.csv',group))
  
  if (selectPerformance & group=='exposure') {
    blinks <- load.DownloadDataframe(informationURLs['blinkdetect'],'blinkdetect_exposure.csv')
    OKparticipants <- blinks$participant[which(blinks$rotated_b == 1 & blinks$performance > 0.65)]
    df <- df[which(df$participant %in% OKparticipants),]
  }
  
  df <- aspligned(df)
  
  # we can select half the localization data... for the exposure group only
  if (group == 'exposure') {
    # print(str(groupdf))
    if (LRpart == 'first') {
      df <- df[which(df$iteration < 3),]
    }
    if (LRpart == 'second') {
      df <- df[which(df$iteration > 2),]
    }
    # print(dim(groupdf))
  } else {
    if (LRpart != 'all') {
      cat('\nWARNING: Partial localization can not be returned for classic/online data.\n\n')
      # the full dataset is returned anyway?
    }
  }
  
  groupdf <- NA
  
  # only get the movement types we need:
  movementnumbers <- c(0,1)
  if (movementtype == 'active') {
    movementnumbers <- c(0)
  }
  if (movementtype == 'passive') {
    movementnumbers <- c(1)
  }
  
  participants <- unique(df$participant)
  
  # remove participants with no data
  removepps <- c()
  for (pp.no in c(1:length(participants))) {
    pp.id <- participants[pp.no]
    for (rotated in c(0,1)) {
      for (passive in movementnumbers) {
        subrows <- which(df$participant == pp.id & df$rotated_b == rotated & df$passive_b == passive)
        # could set this to e.g.: '< 10' to remove participants with very little data?
        if (length(subrows) == 0) {
          removepps <- c(removepps, pp.id)
        }
      }
    }
  }
  if (length(removepps) > 0) {
    participants <- participants[-which(participants %in% unique(removepps))]
    cat('\nremoved these participants:\n')
    print(unique(removepps))
  }
  
  for (pp.no in c(1:length(participants))) {
    
    pp.id <- participants[pp.no]
    
    for (rotated in c(0,1)) {
      
      for (passive in movementnumbers) {
        
        subdf <- df[which(df$participant == pp.id & df$rotated_b == rotated & df$passive_b == passive),]
        
        if (nrow(subdf) == 0) {
          # participant has no data in sub-condition
          next()
        }
        
        locdf <- getLocalizationPoints(subdf, points=points, removeOutliers=TRUE)
        
        if (any(is.na(locdf$taperror_deg))) {
          if (verbose) {
            cat(sprintf('WARNING: NAs in %s, pp:%s, %s/%s\n',group,pp.id,c('aligned','rotated')[rotated+1],c('active','passive')[passive+1]))
          }
        }
        
        if (is.data.frame(groupdf)) {
          groupdf <- rbind(groupdf, locdf)
        } else {
          groupdf <- locdf
        }
        
      }
      
    }
    
  }
  
  # get the dataframe in a useful shape:
  groupdf <- aggregate(taperror_deg ~ group + participant + passive_b + handangle_deg + rotated_b, data=groupdf, FUN=mean, drop=FALSE)
  
  # get rotated - aligned difference, if required:
  if (difference) {
    groupdf <- aggregate(taperror_deg ~ group + participant + passive_b + handangle_deg, data=groupdf, FUN=diff, drop=FALSE)
    groupdf$taperror_deg <- as.numeric(groupdf$taperror_deg)
  }
  
  
  
  # select only one movement type, if required:
  if (movementtype == 'active') {
    groupdf <- groupdf[which(groupdf$passive_b == 0),]
    groupdf <- groupdf[,-which(names(groupdf) == c("passive_b"))]
  }
  if (movementtype == 'passive') {
    groupdf <- groupdf[which(groupdf$passive_b == 1),]
    groupdf <- groupdf[,-which(names(groupdf) == c("passive_b"))]
  }
  
  
  
  return(groupdf)
  
}

aspligned <- function(df) {
  
  participants <- unique(df$participant)
  
  for (pp.no in c(1:length(participants))) {
    
    # fit aligned data with a smooth spline:
    al.idx <- which(df$participant == participants[pp.no] & df$rotated == 0)
    X <- df$handangle_deg[al.idx]
    Y <- df$taperror_deg[al.idx]
    # spar determines smoothness:
    # if set too smooth, some effects disappear
    # print(length(X))
    AL.spl <- smooth.spline(x=X,y=Y,spar=0.90, keep.data=FALSE)
    
    # subtract predicted errors from both aligned and rotated data, using the fitted spline
    pp.idx <- which(df$participant == participants[pp.no])
    predicted_errors <- predict(AL.spl, x=df$handangle_deg[pp.idx])
    df$taperror_deg[pp.idx] <- df$taperror_deg[pp.idx] - predicted_errors$y
    
  }
  
  return(df)
  
}

getLocalizationPoints <- function(subdf, points=c(15,25,35,45,55,65,75), removeOutliers=FALSE) {
  
  X <- subdf$handangle_deg
  Y <- subdf$taperror_deg
  
  # remove data points out of range:
  idx <- which(X > 0 & X < 90)
  X   <- X[idx]
  Y   <- Y[idx]
  
  # see if samples are reasonably close to what a smoothed spline on the rest of the data would predict:
  if (removeOutliers) {
    
    deviations <- c()
    
    for (idx.idx in c(1:length(X))) {
      
      sampleX <- X[idx.idx]
      sampleY <- Y[idx.idx]
      
      # spar determines smoothness:
      # if set too smooth, some effects disappear, but set as high as still seems to work
      spl <- smooth.spline(X[-idx.idx],Y[-idx.idx], spar=.90, keep.data=FALSE) # spar=.90
      
      sampleP <- predict(spl, x=sampleX)
      deviations <- c(deviations, sampleP$y - sampleY)
      
    }
    
    # instead of the mean, we compare each value to its predicted value
    idx <- which(abs(deviations) < (3 * sqrt(sum(deviations^2))))
    
    X <- X[idx]
    Y <- Y[idx]
    
  }
  
  # cubic spline object, based on the data:
  spl <- smooth.spline(X, Y, spar=0.65, keep.data=FALSE) # spar=0.65 # smoothness as low as works
  
  # predict (interpolate) at given coordinates:
  spl.pr <- predict(spl,points)
  # spl.pr$y now has the estimated values of Y (the localization error) at the points of interest
  
  group         <- rep(subdf$group[1],length(points))
  online_b      <- rep(subdf$online_b[1],length(points))
  participant   <- rep(subdf$participant[1],length(points))
  rotated_b     <- rep(subdf$rotated_b[1],length(points))
  passive_b     <- rep(subdf$passive_b[1],length(points))
  handangle_deg <- points
  taperror_deg  <- round(spl.pr$y, digits=5)
  
  # since splines do well at interpolating, but not extrapolating
  # everything outside of the range of the data should be set to NAs:
  Xrange <- range(X)
  
  taperror_deg[which(points < Xrange[1])] <- NA
  taperror_deg[which(points > Xrange[2])] <- NA
  
  locpointdata <- data.frame(group, online_b, participant, rotated_b, passive_b, handangle_deg, taperror_deg)
  locpointdata$participant <- as.character(locpointdata$participant)
  
  return(locpointdata)
  
}

# handle no-cursor data ------

# for every group, this loads the no-cursor reach directions in all relevant tasks,
# and calcuates the reach aftereffects from them
getReachAftereffects <- function(group, part='all', clean=TRUE, difference=TRUE, selectPerformance=TRUE) {
  
  if (group == 'online') {
    group <- 'classic' # same participants, same data, so only stored in one file
  }
  
  # load pre-processed data for the required group (no default)
  raw.df <- load.DownloadDataframe(url=nocursorURLs[group],filename=sprintf('nocursor_%s.csv',group))
  
  if (group=='exposure' & selectPerformance) {
    blinks <- load.DownloadDataframe(informationURLs['blinkdetect'],'blinkdetect_exposure.csv')
    OKparticipants <- blinks$participant[which(blinks$rotated_b == 1 & blinks$performance > 0.65)]
    raw.df <- raw.df[which(raw.df$participant %in% OKparticipants),]
  }
  
  # remove outliers if requested (default: yes)
  if (clean) {
    raw.df <- removeOutliers(raw.df)
  }
  
  # select part of the data if requested (default: use all of it)
  if (part == 'initial') {
    raw.df <- rbind(raw.df[which(raw.df$rotated == 0),], raw.df[which(raw.df$rotated == 1 & raw.df$repetition == 0),])
  }
  if (part == 'remainder') {
    raw.df <- rbind(raw.df[which(raw.df$rotated == 0),], raw.df[which(raw.df$rotated == 1 & raw.df$repetition > 0),])
  }
  
  # average across the task iterations and target repetitions
  avg.df <- aggregate(endpoint_angle ~ participant + rotated + target, data=raw.df, FUN=mean)
  
  # get the difference if requested (default: yes)
  if (difference) {
    avg.df <- aggregate(endpoint_angle ~ participant + target, data=avg.df, FUN=diff)
  }
  
  # return the result:
  return(avg.df)
  
}


removeOutliers <- function(df, stds=2) {
  
  NOKidx <- c()
  OKidx <- c()
  targets <- unique(df$target)
  participants <- unique(df$participant)
  
  for (participant in participants) {
    
    for (rotated in c(0,1)) {
      
      for (target in targets) {
        
        subidx <- which(df$participant == participant & df$target == target & df$rotated == rotated)
        angles <- df$endpoint_angle[subidx]
        # OKidx <- c(OKidx, which(abs(angles - mean(angles)) < (stds * sd(angles))))
        NOKidx <- c(NOKidx, which(abs(angles - mean(angles)) > (stds * sd(angles))))
      }
      
    }
    
  }
  
  Nobs <- nrow(df)
  Nkept <- Nobs - length(NOKidx)
  cat(sprintf('removed %d outliers, kept %0.1f%%\n', Nobs-Nkept, (100 * (Nkept/Nobs))))
  
  # df <- df[OKidx,]
  df$endpoint_angle[NOKidx] <- NA
  
  df <- aggregate(endpoint_angle ~ participant + rotated + repetition + target, data=df, FUN=mean, na.rm=TRUE)
  
  return(df)
  
}

psychometricGeneralization <- function(df, makefigure=FALSE, response='RAE') {
  
  # well... psychometric... Ijust use a cumulative normal distribution function
  
  # does everything use the same dependent variables?
  
  if (makefigure) {
    par(mfrow=c(4,6),mar=c(4,4,3,0.5))
  }
  
  participant <- unique(df$participant)
  peak <- c()
  width <- c()
  
  for (pp.id in participant) {
    
    # let's get the group thing first
    Gdf <- df[-which(df$participant == pp.id),]
    if (response == 'RAE') {
      Gcurve <- aggregate(endpoint_angle ~ target, data=Gdf, FUN=median)
      Gcurve <- Gcurve$endpoint_angle
    }
    if (response == 'loc') {
      Gcurve <- aggregate(taperror_deg ~ handangle_deg, data=Gdf, FUN=median)
      Gcurve <- -1 * Gcurve$taperror_deg
    }
    
    
    # get data for participant:
    Pdf <- df[which(df$participant == pp.id),]
    # print(Pdf)
    
    # get independent variable (target/hand location):
    if (response == 'RAE') {
      IV <- Pdf$target
    }
    if (response == 'loc') {
      IV <- Pdf$handangle_deg
    }
    
    
    
    # get dependent variable (endpoint/tap angle):
    if (response == 'RAE') {
      DV <- Pdf$endpoint_angle
    }
    if (response == 'loc') {
      # Pdf <- Pdf[which(!is.na(Pdf$taperror_deg)),]
      # print(str(Pdf))
      DV <- -1 * Pdf$taperror_deg
    }
    
    # memories of a previous, simpler era:
    # DV <- DV - min(DV)
    # DV <- cumsum(DV)
    # DV <- c(0.01, ((DV/max(DV))*0.97)+0.015, 0.99)
    # unfortunately, it is also wrong
    
    # select non-nan points:
    NN.idx <- which(!is.na(Pdf$taperror_deg))
    if (length(NN.idx) > 0) {
      DV <- DV[NN.idx]
      IV <- IV[NN.idx]
    }
    
    # find an offset to bring the DV on average closer to the group:
    DV <- DV - mean(DV - Gcurve, na.rm=TRUE)
    
    # print(matrix(c(IV,DV), byrow=TRUE, nrow=2))
    
    # if (min(DV) < 0) {
    #   # print(length(which(DV<0)))
    #   DV <- (DV - min(DV))
    # }
    DV <- cumsum(c(0,DV,0))
    DV <- DV/max(DV)
    
    # scale to range
    IV <- c(-45,IV,135)
    if (makefigure) {
      plot(IV,DV,main=sprintf('%s',pp.id),xlab='target angle',ylab='cumulative probability',axes=F,ylim=c(0,1))
    }
    
    # fit pnorm to data:
    # pnorm(X,mean=45,sd=15)
    params <- optim(par=c('mu'=45, 'sigma'=15), fn=pnormMSE, method='Nelder-Mead', Q=IV, Y=DV)
    
    #print(params$par)
    
    peak <- c(peak, params$par['mu'])
    width <- c(width, params$par['sigma'])
    
    if (makefigure) {
      
      lines(IV[2:(length(IV)-1)], c(DV[2],diff(DV[2:(length(DV)-1)]))/(2*max(diff(DV[2:(length(DV)-1)]))), col='green')
      
      X <- seq(-30,120)
      lines(X,pnorm(q=X,mean=params$par['mu'],sd=params$par['sigma']),col='red')
      pp.peak <- params$par['mu']
      lines(c(-30,pp.peak,pp.peak),c(.5,.5,0),lty=2,col='blue')
      
      axis(side=1,at=c(15,45,75))
      axis(side=2,at=c(0,.25,.5,.75,1))
      
      # textlabel <- expression(mu sprintf(': %0.1f\n',pp.peak) sigma sprintf(': %0.1f',params$par['sigma']))
      text(-15,0.90,bquote(mu * ": " * .(round(pp.peak,1)))) 
      text(-15,0.75,bquote(sigma * ": " * .(round(params$par['sigma'],1))))
      # text(x=-15,y=0.90,textlabel)
    }
    
  }
  
  return(data.frame(participant,peak,width))
  
}

pnormMSE <- function(par, Q, Y) {
  
  oldw <- getOption("warn")
  options(warn = -1)
  
  return(mean((pnorm(q=Q, mean=par['mu'], sd=par['sigma']) - Y)^2))
  
  options(warn = oldw)
}

getCDFpeaks <- function(group='classic', response='RAE', movementtype='active', makefigure=FALSE) {
  
  df <- NA
  if (response == 'RAE') {
    df <- getReachAftereffects(group, part='all')
  }
  if (response == 'loc') {
    df <- getPointLocalization(group=group, movementtype=movementtype, verbose=F, difference=TRUE)
  }
  if (!is.data.frame(df)) {
    cat('ERROR: don\'t know what data to load. Fix RESPONSE argument?\n')
    return()
  }
  
  newdf <- psychometricGeneralization(df, makefigure=makefigure, response=response)
  
  return(newdf)
  
}

doLocCDFpeakANOVA <- function(test='training') {
  
  peaks <- NA
  
  if (test == 'training') {
    groups <- c('classic','exposure')
  }
  if (test == 'loctime') {
    groups <- c('classic','online')
  }
  
  for (group in groups) {
    
    for (movementtype in c('active','passive')) {
      
      thisPeakDF <- getCDFpeaks(group=group,response='loc',movementtype=movementtype,makefigure=TRUE)
      thisPeakDF$group <- group
      thisPeakDF$movementtype <- movementtype
      
      if (is.data.frame(peaks)) {
        peaks <- rbind(peaks, thisPeakDF)
      } else {
        peaks <- thisPeakDF
      }
      
    }
    
  }
  
  cols <- c('participant','group','movementtype')
  peaks[cols] <- lapply(peaks[cols], factor)
  
  if (test == 'training') {
    print(ezANOVA(data=peaks, wid=participant, dv=peak, within=movementtype, between=group, type=3))
    print(ezANOVA(data=peaks, wid=participant, dv=width, within=movementtype, between=group, type=3))
  }
  if (test == 'loctime') {
    print(ezANOVA(data=peaks, wid=participant, dv=peak, within=c(movementtype, group), type=3))
    print(ezANOVA(data=peaks, wid=participant, dv=width, within=c(movementtype, group), type=3))
  }
  
}

doRAElocCDFpeakANOVA <- function() {
  
  peaks <- NA
  
  groups <- c('classic','exposure')

  for (group in groups) {
    
    for (response in c('RAE','loc')) {
      
      if (response == 'RAE') {
        thisPeakDF <- getCDFpeaks(group=group,response=response,movementtype='active',makefigure=FALSE)
        thisPeakDF$group <- group
        thisPeakDF$response <- response
        
        if (is.data.frame(peaks)) {
          peaks <- rbind(peaks, thisPeakDF)
        } else {
          peaks <- thisPeakDF
        }
      }
      
      if (response == 'loc') {
        
        for (movementtype in c('active','passive')) {
          
          thisPeakDF <- getCDFpeaks(group=group,response=response,movementtype=movementtype,makefigure=FALSE)
          thisPeakDF$group <- group
          thisPeakDF$response <- sprintf('%s%s',movementtype,response)
          
          if (is.data.frame(peaks)) {
            peaks <- rbind(peaks, thisPeakDF)
          } else {
            peaks <- thisPeakDF
          }
          
        }
        
      }
      
    }
    
  }
  
  peaks <- peaks[-which(peaks$participant == 'gb'),]
  
  cols <- c('participant','group','response')
  peaks[cols] <- lapply(peaks[cols], factor)
  
  peakAOV <- ezANOVA(data=peaks, wid=participant, dv=peak, within=response, between=group, type=3, return_aov=TRUE)
  print(peakAOV[1:3])
  widthAOV <- ezANOVA(data=peaks, wid=participant, dv=width, within=response, between=group, type=3, return_aov=TRUE)
  print(widthAOV[1:3])
  
  print(TukeyHSD(ezANOVA(data=peaks, wid=participant, dv=peak, between=c(response, group), type=3, return_aov=TRUE)$aov))
  print(TukeyHSD(ezANOVA(data=peaks, wid=participant, dv=width, between=c(response, group), type=3, return_aov=TRUE)$aov))
  
  print(aggregate(width ~ response + group, data=peaks, FUN=mean))
  
  print(aggregate(peak ~ response + group, data=peaks, FUN=mean))
  
}