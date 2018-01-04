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
  installRequire.Packages(c('foreach','doParallel'))
  
  cores=detectCores()
  # very friendly for the rest of the system
  usecores <- max(1, ceiling((cores[1] / 2) - 1))
  cl <- makeCluster(usecores)
  registerDoParallel(cl)
  
  x <- as.numeric(names(colMeans(data)))
  
  mus <- foreach (iteration=1:bootstraps, .combine=rbind, .export=c('getGaussianFit','GaussianErrors','parGaussian')) %dopar% {
    
    y <- colMeans(data[sample(row.names(data),size=nrow(data),replace=TRUE),])
    
    as.numeric(getGaussianFit(x,y,mu=47.5,sigma=30,scale=10,offset=4)$par['mu'])
    
  }
  
  stopCluster(cl)
  
  return(quantile(mus,probs=c(.025,.05,.10,.50,.90,.95,.975)))
  
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
