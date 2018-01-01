# first make sure we have all necessary packages installed and loaded:
installRequire.Packages <- function(packages) {
  
  installed.list <- rownames(installed.packages())
  
  for (pkg in packages) {
    
    if (!pkg %in% installed.list) {
      install.packages(pkg,dep=TRUE)
    }
    
    require(pkg, character.only=TRUE)
    
  }
  
}

# not yet sure if we'll use Chi-square, or Satterthwaite estimates of F/p-values:
# installRequire.Packages(c('nlme', 'car', 'lme4', 'lmerTest'))

# localization data can be downloaded from OSF:
groupURLs <- c('exposure'='https://osf.io/47fwu/download', 'classic'='https://osf.io/89t7j/download')

# 
getReachAftereffects <- function(group) {
  
  df <- read.csv(url(groupURLs[group]),stringsAsFactors=FALSE)
  
  df <- aggregate(endpoint_angle ~ participant + rotated + target, data=exposure, FUN=median)
  
  RAE <- aggregate(endpoint_angle ~ participant + target, data=df, FUN=diff)
  
  return(RAE)
  
}