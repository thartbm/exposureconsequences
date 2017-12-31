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
installRequire.Packages(c('nlme', 'car', 'lme4', 'lmerTest'))

# localization data can be downloaded from OSF:

exposure <- read.csv(url('https://osf.io/9qfhp/download'))
classic <- read.csv(url('https://osf.io/upw49/download'))
online <- read.csv(url('https://osf.io/wjcgk/download'))

# no-cursor reach data is not uploaded yet (not sure about the format)

