
source('shared.R')

# first make sure we have all necessary packages installed and loaded:
# installRequire.Packages(c())

# not yet sure if we'll use Chi-square, or Satterthwaite estimates of F/p-values:
# installRequire.Packages(c('nlme', 'car', 'lme4', 'lmerTest'))

required.packages = c('nlme','car')
installRequire.Packages(required.packages)

# localization data can be downloaded from OSF:
groupURLs <- c('exposure'='https://osf.io/9s6au/?action=download', 'classic'='https://osf.io/8hm7f/?action=download')

# for every group, this loads the no-cursor reach directions in all relevant tasks,
# and calcuates the reach aftereffects from them
getReachAftereffects <- function(group, part='all', clean=TRUE) {
  
  raw.df <- load.DownloadDataframe(url=groupURLs[group],filename=sprintf('nocursor_%s.csv',group))
  
  if (clean) {
    clean.df <- removeOutliers(raw.df) 
  }
  
  if (part == 'initial') {
    raw.df <- rbind(raw.df[which(raw.df$rotated == 0),], raw.df[which(raw.df$rotated == 1 & raw.df$repetition == 0),])
  }
  if (part == 'remainder') {
    raw.df <- rbind(raw.df[which(raw.df$rotated == 0),], raw.df[which(raw.df$rotated == 1 & raw.df$repetition > 0),])
  }
  
  avg.df <- aggregate(endpoint_angle ~ participant + rotated + target, data=raw.df, FUN=mean)
  
  RAE <- aggregate(endpoint_angle ~ participant + target, data=avg.df, FUN=diff)
  
  return(RAE)
  
}


removeOutliers <- function(df) {
  
  OKidx <- c()
  targets <- unique(df$target)
  
  for (rotated in c(0,1)) {
    
    for (target in targets) {
      
      subidx <- which(df$target == target & df$rotated == rotated)
      angles <- df$endpoint_angle[subidx]
      OKidx <- c(OKidx, which(abs(angles - mean(angles)) < (2 * sd(angles))))
      
    }
    
  }
  
  Nobs <- nrow(df)
  Nkept <- length(OKidx)
  cat(sprintf('removed %d outliers, kept %0.1f%%\n', Nobs-Nkept, (100 * (Nkept/Nobs))))
  
  df <- df[OKidx,]
  
  df <- aggregate(endpoint_angle ~ participant + rotated + repetition + target, data=df, FUN=mean)
  
  return(df)
  
}


# plot both groups reach aftereffects in one figure
plotReachAftereffects <- function(validate=FALSE) {
  
  points <- c(15,25,35,45,55,65,75)
  
  exposure <- getReachAftereffects('exposure',part='initial')
  classic  <- getReachAftereffects('classic',part='initial')
  exposure_rem <- getReachAftereffects('exposure',part='remainder')
  classic_rem  <- getReachAftereffects('classic',part='remainder')
  
  if (validate) {
    
    validExp <- validateReachAftereffects('exposure')
    validExpParticipants <- validExp$participant[which(validExp$RAE == 1)]
    exposure <- exposure[which(exposure$participant %in% validExpParticipants),]

    validCla <- validateReachAftereffects('classic')
    validClaParticipants <- validCla$participant[which(validCla$RAE == 1)]
    classic <- classic[which(classic$participant %in% validClaParticipants),]
    
  }
  
  exposureAVG <- aggregate(endpoint_angle ~ target, data=exposure, FUN=mean)
  classicAVG <- aggregate(endpoint_angle ~ target, data=classic, FUN=mean)
  
  exposureAVGrem <- aggregate(endpoint_angle ~ target, data=exposure_rem, FUN=mean)
  classicAVGrem <- aggregate(endpoint_angle ~ target, data=classic_rem, FUN=mean)
  
  exposureCI <- matrix(unlist(by(exposure$endpoint_angle, INDICES=c(exposure$target), FUN=t.interval)),nrow=2)
  classicCI <- matrix(unlist(by(classic$endpoint_angle, INDICES=c(classic$target), FUN=t.interval)),nrow=2)
  
  X <- c(points, rev(points))
  expY <- c(exposureCI[1,],rev(exposureCI[2,]))
  claY <- c(classicCI[1,],rev(classicCI[2,]))
  
  plot(-1000,-1000, main='reach aftereffects', xlab='target angle [deg]', ylab='reach endpoint deviation [deg]', xlim=c(10,80), ylim=c(0,15), axes=F)
  
  polygon(X,expY,border=NA,col=rgb(1,0,0,.2))
  polygon(X,claY,border=NA,col=rgb(0,0,1,.2))
  
  lines(points,exposureAVG$endpoint_angle,col=rgb(1,0,0))
  lines(points,classicAVG$endpoint_angle,col=rgb(0,0,1))
  lines(points,exposureAVGrem$endpoint_angle,col=rgb(1,0,0),lty=2)
  lines(points,classicAVGrem$endpoint_angle,col=rgb(0,0,1),lty=2)
  
  axis(1,at=points)
  axis(2,at=c(0,5,10,15))
  
  legend(10,15,c('exposure','classic'),col=c(rgb(1,0,0),rgb(0,0,1)),lty=c(1,1),bty='n')
  
}
# comment
getPeakConfidenceInterval <- function(group,validate=FALSE,part='initial') {
  
  RAE <- getReachAftereffects(group, part=part)
  
  if (validate) {
    validParticipants <- validateReachAftereffects(group, method='45')
    validParticipants <- validParticipants$participant[which(validParticipants$RAE == 1)]
    RAE <- RAE[which(RAE$participant %in% validParticipants),]
  }
  
  RAE <- xtabs(endpoint_angle~.,RAE)
  
  bootstrapGaussianPeak(data=RAE,bootstraps=1000,mu=47.5,sigma=30,scale=10,offset=4)
  
}

validateReachAftereffects <- function(group, method='45') {
  
  # method can be: '45', 'all_aov', 'all_ttest', 'all_any5deg'
  
  raw.df <- load.DownloadDataframe(url=groupURLs[group],filename=sprintf('%s_nocursor.csv',group))
  participant <- unique(raw.df$participant)
  
  if (method == '45') {
    aligned <- raw.df[which(raw.df$rotated == 0 & raw.df$target == 45),c('participant','repetition','endpoint_angle')]
    rotated <- raw.df[which(raw.df$rotated == 1 & raw.df$repetition == 0 & raw.df$target == 45),c('participant','endpoint_angle')]
  } else {
    aligned <- aggregate(endpoint_angle ~ target + participant, data=raw.df[which(raw.df$rotated == 0),c('participant','repetition','target','endpoint_angle')], FUN=mean)[,c(2,1,3)]
    rotated <- raw.df[which(raw.df$rotated == 1 & raw.df$repetition == 0),c('participant','target','endpoint_angle')]
  }
  
  if (method == 'all_aov') {
    installRequire.Packages(c('afex'))
  }
  
  RAE <- c()
  
  for (ppno in participant) {
    
    ppAL <- aligned$endpoint_angle[aligned$participant == ppno]
    ppRO <- rotated$endpoint_angle[rotated$participant == ppno]
    
    if (method == '45') {
      pptt <- t.test(ppAL, mu=ppRO, alternative='less')
      if (pptt$p.value < .05) {
        RAE <- c(RAE, 1)
      } else {
        RAE <- c(RAE, 0)
      }
    }
    
    if (method == 'all_aov') {
      session <- c(rep(0,7),rep(1,7))
      target <- c(ppAL$target, ppRO$target)
      endpoint_angle <- c(ppAL$endpoint_angle, ppRO$endpoint_angle)

      ppDF <- data.frame(session,target,endpoint_angle)
      ppDF$session <- factor(ppDF$session)
      ppDF$target <- factor(ppDF$target)
      
      ppEZ <- aov_ez(id='target', dv='endpoint_angle', data=ppDF, within='session')
      if (ppEZ[[1]][1,6] < .05) {
        RAE <- c(RAE, 1)
      } else {
        RAE <- c(RAE, 0)
      }
    }
    
    if (method == 'all_ttest') {
      pptt <- t.test(x=ppAL$endpoint_angle, y=ppRO$endpoint_angle, paired=TRUE, alternative='l')
      if (pptt$p.value < .05) {
        RAE <- c(RAE, 1)
      } else {
        RAE <- c(RAE, 0)
      }
    }
    
    if (method == 'all_any5deg') {
      rae <- ppRO$endpoint_angle - ppAL$endpoint_angle
      if (any(rae > 5)) {
        RAE <- c(RAE, 1)
      } else {
        RAE <- c(RAE, 0)
      }      
    }
    
  }
  
  return(data.frame(participant, RAE))
  
}

plotReachAftereffectDistributions <- function() {
  
  # if exposure training does not engage typical adaptation,
  # but rather model-free learning, previously successful movements
  # would be used, which can be one of two:
  # 1) movements such as those with aligned/normal feedback
  # 2) the movements done during the rotated/exposure training
  #
  # we ignore participants, and plot a kernel-density estimation
  # first of the rotated session only
  
  exposure <- getReachAftereffects('exposure')
  classic  <- getReachAftereffects('classic')
  
  
  
}