
source('shared.R')

# first make sure we have all necessary packages installed and loaded:
# installRequire.Packages(c())

# not yet sure if we'll use Chi-square, or Satterthwaite estimates of F/p-values:
# installRequire.Packages(c('nlme', 'car', 'lme4', 'lmerTest'))

# localization data can be downloaded from OSF:
groupURLs <- c('exposure'='https://osf.io/47fwu/download', 'classic'='https://osf.io/89t7j/download')

# for every group, this loads the no-cursor reach directions in all relevant tasks,
# and calcuates the reach aftereffects from them
getReachAftereffects <- function(group) {
  
  raw.df <- read.csv(url(groupURLs[group]),stringsAsFactors=FALSE)
  
  avg.df <- aggregate(endpoint_angle ~ participant + rotated + target, data=raw.df, FUN=mean)
  
  RAE <- aggregate(endpoint_angle ~ participant + target, data=avg.df, FUN=diff)
  
  return(RAE)
  
}

# plot both groups reach aftereffects in one figure
plotReachAftereffects <- function() {
  
  points <- c(15,25,35,45,55,65,75)
  
  exposure <- getReachAftereffects('exposure')
  classic  <- getReachAftereffects('classic')
  
  exposureAVG <- aggregate(endpoint_angle ~ target, data=exposure, FUN=mean)
  classicAVG <- aggregate(endpoint_angle ~ target, data=classic, FUN=mean)
  
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
  
  axis(1,at=points)
  axis(2,at=c(0,5,10,15))
  
  legend(10,15,c('exposure','classic'),col=c(rgb(1,0,0),rgb(0,0,1)),bty='n')
  
}