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
  
  exposure <- getReachAftereffects('exposure')
  classic  <- getReachAftereffects('classic')
  
  
  
}