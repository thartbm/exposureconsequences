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

bootstrapGaussianPeak <- function(data,bootstraps=1000,mu=47.5,sigma=30,scale=10,offset=4) {
  
  # parallel for-loop?
  # kernels <- max(1, ceil((total kernels / 2) - 1)) # very friendly for the rest of the system
  
  mus <- c()
  
  x <- as.numeric(names(colMeans(data)))
  
  for (iteration in c(1:bootstraps)) {
    
    y <- colMeans(data[sample(row.names(data),size=nrow(data),replace=TRUE),])
    
    mus <- c(mus,as.numeric(getGaussianFit(x,y,mu=47.5,sigma=30,scale=10,offset=4)$par['mu']))
    
  }
  
  return(quantile(mus,probs=c(.025,.50,.975)))
  
}

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
