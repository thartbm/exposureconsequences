# some functions shared by the analyses of the two types of data

installRequire.Packages <- function(packages) {
  
  installed.list <- rownames(installed.packages())
  
  for (pkg in packages) {
    
    if (!pkg %in% installed.list) {
      install.packages(pkg,dep=TRUE)
    }
    
    require(pkg, character.only=TRUE)
    
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