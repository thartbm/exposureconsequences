# some functions shared by the analyses of the two types of data

# localization data can be downloaded from OSF:
nocursorURLs <- c('exposure'='https://osf.io/9s6au/?action=download', 'classic'='https://osf.io/8hm7f/?action=download')

# no-cursor data can be downloaded from OSF:
localizationURLs <- c('exposure'='https://osf.io/9f6gu/?action=download', 'classic'='https://osf.io/upw49/?action=download', 'online'='https://osf.io/wjcgk/download')

# we'll control the color in figures centrally from here:
colorset <- list()

colorset[['expActS']] <- '#005de4ff' # blue
colorset[['expActT']] <- '#005de42f'
colorset[['expPasS']] <- '#2ab2f2ff' # lighter blue
colorset[['expPasT']] <- '#2ab2f22f'
colorset[['claActS']] <- '#e51636ff' # "York red"
colorset[['claActT']] <- '#e516362f'
colorset[['claPasS']] <- '#ff8200ff' # orange
colorset[['claPasT']] <- '#ff82002f'

colorset[['extra1S']] <- '#c400c4ff' # purple
colorset[['extra1T']] <- '#c400c42f'


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
    cat('\ncar and nlme are required to reproduce the analyses\n')
    cat('\nforeach and doParallel speed up a few functions\n')
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

getPointLocalization <- function(group, difference=TRUE, points=c(15,25,35,45,55,65,75), movementtype='both', verbose=TRUE, LRpart='all') {
  
  df <- load.DownloadDataframe(url=localizationURLs[group],filename=sprintf('localization_%s.csv',group))
  
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
    if (LRpart != 'both') {
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
  X <- X[idx]
  Y <- Y[idx]
  
  # see if samples are reasonably close to what a smoothed spline on the rest of the data would predict:
  if (removeOutliers) {
    
    deviations <- c()
    
    for (idx.idx in c(1:length(X))) {
      
      sampleX <- X[idx.idx]
      sampleY <- Y[idx.idx]
      
      # spar determines smoothness:
      # if set too smooth, some effects disappear
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
  spl <- smooth.spline(X, Y, spar=0.65, keep.data=FALSE) # spar=0.65
  
  # predict (interpolate) at given coordinates:
  spl.pr <- predict(spl,points)
  # spl.pr$y now has the estimated values of Y (the localization error) at the points of interest
  
  group <- rep(subdf$group[1],length(points))
  online_b <- rep(subdf$online_b[1],length(points))
  participant <- rep(subdf$participant[1],length(points))
  rotated_b <- rep(subdf$rotated_b[1],length(points))
  passive_b <- rep(subdf$passive_b[1],length(points))
  handangle_deg <- points
  taperror_deg <- round(spl.pr$y, digits=5)
  
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
getReachAftereffects <- function(group, part='all', clean=TRUE, difference=TRUE) {
  
  if (group == 'online') {
    group <- 'classic'
  }
  
  # load pre-processed data for the required group (no default)
  raw.df <- load.DownloadDataframe(url=nocursorURLs[group],filename=sprintf('nocursor_%s.csv',group))
  
  # remove outliers if requested (default: yes)
  if (clean) {
    clean.df <- removeOutliers(raw.df) 
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
  
  # return te result:
  return(avg.df)
  
}



removeOutliers <- function(df, stds=2) {
  
  OKidx <- c()
  targets <- unique(df$target)
  participants <- unique(df$participant)
  
  for (participant in participants) {
    
    for (rotated in c(0,1)) {
      
      for (target in targets) {
        
        subidx <- which(df$participant == participant & df$target == target & df$rotated == rotated)
        angles <- df$endpoint_angle[subidx]
        OKidx <- c(OKidx, which(abs(angles - mean(angles)) < (stds * sd(angles))))
        
      }
      
    }
    
  }
  
  Nobs <- nrow(df)
  Nkept <- length(OKidx)
  cat(sprintf('removed %d outliers, kept %0.1f%%\n', Nobs-Nkept, (100 * (Nkept/Nobs))))
  
  df <- df[OKidx,]
  
  df <- aggregate(endpoint_angle ~ participant + rotated + repetition + target, data=df, FUN=mean)
  
  return(df)
  
}

