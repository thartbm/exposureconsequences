
source('shared.R')

correlateNoCursorsLocalization <- function(NCpart='all', generateSVG=FALSE, selectPerformance=TRUE) {
  
  # first we do all responses, plan to split it by target later
  # get all data, for classic and exposure, we need localization and reach aftereffects
  
  # LOCALIZATION SHIFT for exposure
  local.exp <- getPointLocalization(group='exposure', difference=FALSE, verbose=FALSE, selectPerformance=selectPerformance)
  # active
  loc.exp.act.ro <- aggregate(taperror_deg ~ participant + handangle_deg, data=local.exp[local.exp$rotated_b == 1 & local.exp$passive_b == 0,], FUN=mean, drop=FALSE)
  loc.exp.act.al <- aggregate(taperror_deg ~ participant + handangle_deg, data=local.exp[local.exp$rotated_b == 0 & local.exp$passive_b == 0,], FUN=mean, drop=FALSE)
  local.exp.act <- loc.exp.act.ro
  local.exp.act$taperror_deg <- local.exp.act$taperror_deg - loc.exp.act.al$taperror_deg
  # passive
  loc.exp.pas.ro <- aggregate(taperror_deg ~ participant + handangle_deg, data=local.exp[local.exp$rotated_b == 1 & local.exp$passive_b == 1,], FUN=mean, drop=FALSE)
  loc.exp.pas.al <- aggregate(taperror_deg ~ participant + handangle_deg, data=local.exp[local.exp$rotated_b == 0 & local.exp$passive_b == 1,], FUN=mean, drop=FALSE)
  local.exp.pas <- loc.exp.pas.ro
  local.exp.pas$taperror_deg <- local.exp.pas$taperror_deg - loc.exp.pas.al$taperror_deg
  
  # REACH AFTEREFFECTS for exposure
  nocur.exp <- getReachAftereffects(group='exposure', part=NCpart, clean=TRUE, selectPerformance=selectPerformance) 
  names(nocur.exp)[names(nocur.exp) == 'endpoint_angle'] <- 'RAE'
  # combine data frames for exposure
  exposure <- nocur.exp[which(nocur.exp$participant %in% local.exp.act$participant),]
  exposure$active_localization_shift <- local.exp.act$taperror_deg
  exposure$passive_localization_shift <- local.exp.pas$taperror_deg
  
  # LOCALIZATION SHIFT for classic
  local.cla <- getPointLocalization(group='classic', difference=FALSE, verbose=FALSE, selectPerformance=selectPerformance)
  # active
  loc.cla.act.ro <- aggregate(taperror_deg ~ participant + handangle_deg, data=local.cla[local.cla$rotated_b == 1 & local.cla$passive_b == 0,], FUN=mean, drop=FALSE)
  loc.cla.act.al <- aggregate(taperror_deg ~ participant + handangle_deg, data=local.cla[local.cla$rotated_b == 0 & local.cla$passive_b == 0,], FUN=mean, drop=FALSE)
  local.cla.act <- loc.cla.act.ro
  local.cla.act$taperror_deg <- local.cla.act$taperror_deg - loc.cla.act.al$taperror_deg
  # passive
  loc.cla.pas.ro <- aggregate(taperror_deg ~ participant + handangle_deg, data=local.cla[local.cla$rotated_b == 1 & local.cla$passive_b == 1,], FUN=mean, drop=FALSE)
  loc.cla.pas.al <- aggregate(taperror_deg ~ participant + handangle_deg, data=local.cla[local.cla$rotated_b == 0 & local.cla$passive_b == 1,], FUN=mean, drop=FALSE)
  local.cla.pas <- loc.cla.pas.ro
  local.cla.pas$taperror_deg <- local.cla.pas$taperror_deg - loc.cla.pas.al$taperror_deg
  
  # REACH AFTEREFFECTS for classic
  nocur.cla <- getReachAftereffects(group='classic', part=NCpart, clean=TRUE, selectPerformance=selectPerformance) 
  names(nocur.cla)[names(nocur.cla) == 'endpoint_angle'] <- 'RAE'
  # combine the dataframes for classic:
  classic <- nocur.cla[which(nocur.cla$participant %in% local.cla$participant),]
  classic$active_localization_shift <- local.cla.act$taperror_deg
  classic$passive_localization_shift <- local.cla.pas$taperror_deg
  
  # GOT DATA, NOW ANALYZE... oh wait
  
  # the correlations are polluted by the generalization curve... let's remove that!
  
  # get the average across a set of targets (the flat area? what do we need/want?)
  classic.tmp <- classic[which(classic$target %in% c(35,45,55,65)),]
  classic <- aggregate(RAE ~ participant, data=classic.tmp, FUN=mean, na.rm=TRUE, drop=FALSE)
  classic$active_localization_shift <- aggregate(active_localization_shift ~ participant, data=classic.tmp, FUN=mean, na.rm=TRUE, drop=FALSE)$active_localization_shift
  classic$passive_localization_shift <- aggregate(passive_localization_shift ~ participant, data=classic.tmp, FUN=mean, na.rm=TRUE, drop=FALSE)$passive_localization_shift
  # and for exposure too
  exposure.tmp <- exposure[which(exposure$target %in% c(35,45,55,65)),]
  exposure <- aggregate(RAE ~ participant, data=exposure.tmp, FUN=mean, na.rm=TRUE, drop=FALSE)
  exposure$active_localization_shift <- aggregate(active_localization_shift ~ participant, data=exposure.tmp, FUN=mean, na.rm=TRUE, drop=FALSE)$active_localization_shift
  exposure$passive_localization_shift <- aggregate(passive_localization_shift ~ participant, data=exposure.tmp, FUN=mean, na.rm=TRUE, drop=FALSE)$passive_localization_shift
  
  
  # # NOW do the tests:
  xlimits <- rev(c(-20.5,10.5))
  ylimits <- c(-10.5,20.5)
  xticks <- seq(from=-20,to=10,by=10)
  yticks <- rev(seq(from=-10,to=20,by=10))
  
  ylab <- 'reach aftereffect [°]'
  xlab <- 'localization shift [°]'
  
  if (generateSVG) {
    installed.list <- rownames(installed.packages())
    if ('svglite' %in% installed.list) {
      library('svglite')
      svglite(file='Fig6.svg', width=7.5, height=2.5, system_fonts=list(sans='Arial', mono='Times New Roman'))
    } else {
      generateSVG=FALSE
    }
  }
  
  par(mfrow=c(1,3))
  
  dataframes = list()
  dataframes[['exposure']] <- exposure
  dataframes[['classic']]  <- classic
  
  output <- list()
  
  groups <- c('exposure','classic')
  
  for (groupno in c(1:length(groups))) {
    
    group <- groups[groupno]
    
    plot(-1000, -1000, main=group, xlab=xlab, ylab=ylab, asp=1, axes=FALSE, xlim=xlimits, ylim=ylimits)
    
    mtext(toupper(letters[groupno]), side=3, outer=TRUE, at=c((groupno-1)/3,1), line=-1, adj=0, padj=1)
    
    lines(c(0,0),range(yticks),col=rgb(.5,.5,.5,.5))
    lines(range(xticks),c(0,0),col=rgb(.5,.5,.5,.5))
    lines(rev(range(xticks)),range(yticks),col=rgb(.5,.5,.5,.5),lty=3)
    
    legendlinecolors <- c()
    
    for (movement in c('active','passive')) {
      df <- dataframes[[group]]
      colname <- sprintf('%s_localization_shift',movement)
      X <- df[,colname]
      Y <- df$RAE
      
      colorID <- sprintf('%s%s%s',substring(group,1,3),toupper(substring(movement,1,1)),substring(movement,2,3))
      Solid       <- colorset[[sprintf('%sS',colorID)]]
      Transparent <- colorset[[sprintf('%sT',colorID)]]
      
      legendlinecolors <- c(legendlinecolors,Solid)
      
      points(X,Y,col=Solid,pch=1,cex=1.2)
      
      thisoutput <- list()
      thisoutput[['cortest']] <- cor.test(X,Y)
      thisoutput[['linearm']] <- summary(lm(Y ~ X))
      output[[toupper(sprintf('%s %s\n',group,movement))]] <- thisoutput
      
      colors = c(Transparent,Solid)
      
      plotCorrelationWithCI(X=X,Y=Y,colors=colors)
      
    }
    
    legend(-3,-3,c('active','passive'),col=legendlinecolors,lwd=c(1.5,1.5),bty='n')
    
    if (length(xticks)) {
      axis(1,at=xticks)
    }
    if (length(yticks)) {
      axis(2,at=yticks)
    }
    
  }
  
  if (generateSVG) {
    dev.off()
  }
  
  return(output)
  
}


multipleRegressionLocalization <- function(group,NCpart='all',expart='all',LRpart='all', selectPerformance=TRUE) {
  
  # first we do all responses, plan to split it by target later
  # get all data, for classic and exposure, we need localization and reach aftereffects
  
  # LOCALIZATION SHIFT
  local <- getPointLocalization(group=group, difference=FALSE, verbose=FALSE, LRpart=LRpart, selectPerformance=selectPerformance)
  # active
  loc.act.ro <- aggregate(taperror_deg ~ participant + handangle_deg, data=local[local$rotated_b == 1 & local$passive_b == 0,], FUN=mean, drop=FALSE)
  loc.act.al <- aggregate(taperror_deg ~ participant + handangle_deg, data=local[local$rotated_b == 0 & local$passive_b == 0,], FUN=mean, drop=FALSE)
  local.act <- loc.act.ro
  local.act$taperror_deg <- local.act$taperror_deg - loc.act.al$taperror_deg
  # passive
  loc.pas.ro <- aggregate(taperror_deg ~ participant + handangle_deg, data=local[local$rotated_b == 1 & local$passive_b == 1,], FUN=mean, drop=FALSE)
  loc.pas.al <- aggregate(taperror_deg ~ participant + handangle_deg, data=local[local$rotated_b == 0 & local$passive_b == 1,], FUN=mean, drop=FALSE)
  local.pas <- loc.pas.ro
  local.pas$taperror_deg <- local.pas$taperror_deg - loc.pas.al$taperror_deg
  
  # REACH AFTEREFFECTS
  nocur <- getReachAftereffects(group=group, part=NCpart, clean=TRUE, selectPerformance=selectPerformance) 
  names(nocur)[names(nocur) == 'endpoint_angle'] <- 'RAE'
  # combine data frames
  df <- nocur[which(nocur$participant %in% local.act$participant),]
  df$active <- local.act$taperror_deg
  df$passive <- local.pas$taperror_deg
  

  # the correlations are polluted by the generalization curve... let's remove that!
  
  # and for exposure too
  df.tmp <- df[which(df$target %in% c(35,45,55,65)),]
  df <- aggregate(RAE ~ participant, data=df.tmp, FUN=mean, na.rm=TRUE, drop=FALSE)
  df$active <- aggregate(active ~ participant, data=df.tmp, FUN=mean, na.rm=TRUE, drop=FALSE)$active
  df$passive <- aggregate(passive ~ participant, data=df.tmp, FUN=mean, na.rm=TRUE, drop=FALSE)$passive
  
  cat(sprintf('\n%s\n\n',toupper(group)))
  
  # str(df)
  
  # THIS CAUSES AN ERROR NOW: data argument is of the wrong type...
  # because STEP is imported from lmerTest, but we want the basic one from 'stats'
  fm <- stats::step(lm(RAE ~ active + passive, data=df))
  
  print(summary(fm))
  
}


plotCorrelationWithCI <- function(X,Y,colors=c('black','red'),main='main title',xlab='x',ylab='y',axes=FALSE,xticks=c(),yticks=c(),xlim=NULL,ylim=NULL) {
  
  # fit regression model
  this.lm <- lm(Y ~ X)
  
  # now show the confidence interval
  pointlocs <- seq(min(X),max(X),.1)
  
  # lines(
  #   x   = pointlocs,
  #   y   = predict( this.lm, 
  #                  newdata=data.frame(X=pointlocs), 
  #                  interval = "confidence" )[ , "upr" ],
  #   col = colors[2],
  #   lty = 2)
  # 
  # lines(
  #   x   = pointlocs,
  #   y   = predict( this.lm, 
  #                  newdata=data.frame(X=pointlocs), 
  #                  interval = "confidence" )[ , "lwr" ],
  #   col = colors[2],
  #   lty = 2)
  
  y1 = predict( this.lm, newdata=data.frame(X=pointlocs), interval = "confidence" )[ , "upr" ]
  y2 = predict( this.lm, newdata=data.frame(X=pointlocs), interval = "confidence" )[ , "lwr" ]
  
  polygon(c(pointlocs,rev(pointlocs)),c(y1,rev(y2)), col=colors[1], border=NA)
  
  # and use that to show a line:
  lines(range(X), predict(this.lm, newdata=data.frame(X=range(X))), col=colors[2], lwd=1.5)
  
}

# Generalization -----

plotGeneralization <- function(generateSVG=FALSE, selectPerformance=TRUE, remove15=FALSE) {
  
  if (generateSVG) {
    installed.list <- rownames(installed.packages())
    if ('svglite' %in% installed.list) {
      library('svglite')
      svglite(file='Fig5.svg', width=7.5, height=3, system_fonts=list(sans='Arial', mono='Times New Roman'))
    } else {
      generateSVG=FALSE
    }
  }
  
  
  par(mfrow=c(1,3))
  
  # get the data to plot:
  exp <- getPointLocalization('exposure', difference=TRUE, verbose=FALSE, selectPerformance=selectPerformance)
  cla <- getPointLocalization('classic', difference=TRUE, verbose=FALSE, selectPerformance=FALSE)
  
  # get the averages for the line plots:
  exp.avg <- aggregate(taperror_deg ~ passive_b + handangle_deg, data=exp, FUN=mean)
  exp.avg.act <- exp.avg[which(exp.avg$passive_b == 0),]
  exp.avg.pas <- exp.avg[which(exp.avg$passive_b == 1),]
  cla.avg <- aggregate(taperror_deg ~ passive_b + handangle_deg, data=cla, FUN=mean)
  cla.avg.act <- cla.avg[which(cla.avg$passive_b == 0),]
  cla.avg.pas <- cla.avg[which(cla.avg$passive_b == 1),]
  
  points=c(15,25,35,45,55,65,75)
  exp <- getPointLocalization(group='exposure', difference=TRUE, points=points, movementtype='both', LRpart='all', verbose=FALSE, selectPerformance=selectPerformance)
  cla <- getPointLocalization(group='classic', difference=TRUE, points=points, movementtype='both', LRpart='all', verbose=FALSE, selectPerformance=selectPerformance)
  
  exp.act <- exp[which(exp$passive_b == 0 & is.finite(exp$taperror_deg)),]
  cla.act <- cla[which(cla$passive_b == 0 & is.finite(cla$taperror_deg)),]
  
  # cla.avg.act <- aggregate(taperror_deg ~ handangle_deg, data=cla.act, FUN=mean) 
  # exp.avg.act <- aggregate(taperror_deg ~ handangle_deg, data=exp.act, FUN=mean) 
  
  # if (remove15) {
  #   exp.act <- exp.act[which(exp.act$handangle_deg > 15),]
  #   cla.act <- cla.act[which(cla.act$handangle_deg > 15),]
  # }
  
  # fitting on all data:
  exp.fit <- getGaussianFit(x=exp.act$handangle_deg,exp.act$taperror_deg,mu=50,sigma=10,scale=-75,offset=-4)
  cla.fit <- getGaussianFit(x=cla.act$handangle_deg,cla.act$taperror_deg,mu=50,sigma=10,scale=-75,offset=-4)
  
  # get confidence intervals for the peak of the generalization curve for localization shifts:
  cla.locshift <- getPeakLocConfidenceInterval(group='classic',
                                               CIs=c(.95), 
                                               movementtype='active', 
                                               LRpart='all',
                                               selectPerformance=FALSE,
                                               remove15=remove15)
  exp.locshift <- getPeakLocConfidenceInterval(group='exposure',
                                               CIs=c(.95), 
                                               movementtype='active', 
                                               LRpart='all', 
                                               selectPerformance=selectPerformance,
                                               remove15=remove15)
  
  plot(-1000,-1000, main='localization shifts', xlab='hand angle [°]', ylab='localization shift [°]', xlim=c(10,80), ylim=c(0,-15), axes=F)
  
  mtext('A', side=3, outer=TRUE, at=c(0,1), line=-1, adj=0, padj=1)
  
  # plot the data, faintly
  
  # lines(points[1:2],cla.avg.act$taperror_deg[1:2],col=colorset[['claActT']],lty=2,lwd=1.5)
  # lines(points[2:7],cla.avg.act$taperror_deg[2:7],col=colorset[['claActT']],lty=1,lwd=1.5)
  # lines(points[1:2],exp.avg.act$taperror_deg[1:2],col=colorset[['expActT']],lty=2,lwd=1.5)
  # lines(points[2:7],exp.avg.act$taperror_deg[2:7],col=colorset[['expActT']],lty=1,lwd=1.5)
  
  lines(points,cla.avg.act$taperror_deg,col=colorset[['claActT']],lty=1,lwd=1.5)
  lines(points,exp.avg.act$taperror_deg,col=colorset[['expActT']],lty=1,lwd=1.5)
  
  # plot fitted Gaussian functions to all data:
  
  X <- seq(15,75)
  cla.Y.fit <- cla.fit$par['scale']*parGaussian(cla.fit$par,X)
  cla.Y.fit <- cla.Y.fit + cla.fit$par['offset']
  exp.Y.fit <- exp.fit$par['scale']*parGaussian(exp.fit$par,X)
  exp.Y.fit <- exp.Y.fit + exp.fit$par['offset']
  
  lines(X,cla.Y.fit,col=colorset[['claActS']],lty=1,lw=1.5)
  lines(X,exp.Y.fit,col=colorset[['expActS']],lty=1,lw=1.5)
  
  cla.idx <- which.min(cla.Y.fit)
  exp.idx <- which.min(exp.Y.fit)
  
  # connect peaks of group fits to CIs:
  
  arrows(X[cla.idx],cla.Y.fit[cla.idx],X[cla.idx],-2.5,col=colorset[['claActS']],lwd=1.5,length=.05)
  arrows(X[exp.idx],exp.Y.fit[exp.idx],X[exp.idx],-2.5,col=colorset[['expActS']],lwd=1.5,length=.05)
  
  # indicate feedback and hand position during training:
  
  arrows(45,-2.5,45,-1,col='black',lw=1.5,length=0.05)
  arrows(75,0,75,-1.5,col='black',lw=1.5,length=0.05)
  
  # plot the bootstrap peaks of the generalization functions
  
  polygon(cla.locshift$value[c(1,3,3,1)],c(0,0,-1,-1),border=NA,col=colorset[['claActT']])
  polygon(exp.locshift$value[c(1,3,3,1)],c(-1.5,-1.5,-2.5,-2.5),border=NA,col=colorset[['expActT']])
  
  lines(cla.locshift$value[c(2,2)],c(0,-2.5),col=colorset[['claActS']],lty=1,lw=1.5)
  lines(exp.locshift$value[c(2,2)],c(0,-2.5),col=colorset[['expActS']],lty=1,lw=1.5)
  
  # add legend
  legend(10,-15,c('exposure','classic'),col=c(colorset[['expActS']],colorset[['claActS']]),lty=c(1,1),lwd=c(1.5,1.5),bty='n')
  
  # add tick marks:
  axis(1,at=points)
  axis(2,at=c(0,-5,-10,-15))
  
  
  
  # B      Reach Aftereffects
  
  
  classic  <- getReachAftereffects('classic',part='all', selectPerformance=selectPerformance)
  exposure <- getReachAftereffects('exposure',part='all', selectPerformance=selectPerformance)
  
  # exposureAVGini <- aggregate(endpoint_angle ~ target, data=exposure_ini, FUN=mean) # error?
  # exposureAVGrem <- aggregate(endpoint_angle ~ target, data=exposure_rem, FUN=mean)
  classicAVG <- aggregate(endpoint_angle ~ target, data=classic, FUN=mean)
  exposureAVG <- aggregate(endpoint_angle ~ target, data=exposure, FUN=mean)
  
  
  points=c(15,25,35,45,55,65,75)
  
  classic  <- getReachAftereffects('classic',part='all', selectPerformance=selectPerformance)
  exposure <- getReachAftereffects('exposure',part='all', selectPerformance=selectPerformance)
  
  # fitting group data:
  exp.fit <- getGaussianFit(x=exposure$target,exposure$endpoint_angle,mu=50,sigma=30,scale=50,offset=4)
  cla.fit <- getGaussianFit(x=classic$target,classic$endpoint_angle,mu=50,sigma=30,scale=50,offset=4)
  
  # get confidence intervals for the peak of the generalization curve for localization shifts:
  cla.RAEshift <- getPeakConfidenceInterval('classic', part='all', CIs=c(.95), selectPerformance=selectPerformance)
  exp.RAEshift <- getPeakConfidenceInterval('exposure', part='all', CIs=c(.95), selectPerformance=selectPerformance)
  
  plot(-1000,-1000, main='reach aftereffects', xlab='target angle [°]', ylab='reach endpoint deviation [°]', xlim=c(10,80), ylim=c(0,15), axes=F)
  
  mtext('B', side=3, outer=TRUE, at=c(1/3,1), line=-1, adj=0, padj=1)
  
  # plot the data, faintly
  
  lines(points,classicAVG$endpoint_angle,col=colorset[['claActT']],lwd=1.5)
  lines(points,exposureAVG$endpoint_angle,col=colorset[['expActT']],lwd=1.5)
  
  # plot fitted Gaussian functions to all data:
  
  X <- seq(15,75)
  cla.Y.fit <- cla.fit$par['scale']*parGaussian(cla.fit$par,X)
  cla.Y.fit <- cla.Y.fit + cla.fit$par['offset']
  exp.Y.fit <- exp.fit$par['scale']*parGaussian(exp.fit$par,X)
  exp.Y.fit <- exp.Y.fit + exp.fit$par['offset']
  
  lines(X,cla.Y.fit,col=colorset[['claActS']],lty=1,lw=1.5)
  lines(X,exp.Y.fit,col=colorset[['expActS']],lty=1,lw=1.5)
  
  cla.idx <- which.max(cla.Y.fit)
  exp.idx <- which.max(exp.Y.fit)
  
  # connect peaks of group fits to CIs:
  
  arrows(X[cla.idx],cla.Y.fit[cla.idx],X[cla.idx],2.5,col=colorset[['claActS']],lwd=1.5,length=.05)
  arrows(X[exp.idx],exp.Y.fit[exp.idx],X[exp.idx],2.5,col=colorset[['expActS']],lwd=1.5,length=.05)
  
  # indicate feedback and hand position during training:
  
  arrows(45,2.5,45,1,col='black',lw=1.5,length=0.05)
  arrows(75,0,75,1.5,col='black',lw=1.5,length=0.05)
  
  # plot the bootstrap peaks of the generalization functions
  
  polygon(cla.RAEshift$value[c(1,3,3,1)],c(0,0,1,1),border=NA,col=colorset[['claActT']])
  polygon(exp.RAEshift$value[c(1,3,3,1)],c(1.5,1.5,2.5,2.5),border=NA,col=colorset[['expActT']])
  
  lines(cla.RAEshift$value[c(2,2)],c(0,2.5),col=colorset[['claActS']],lty=1,lw=1.5)
  lines(exp.RAEshift$value[c(2,2)],c(0,2.5),col=colorset[['expActS']],lty=1,lw=1.5)
  
  # add tick marks:
  axis(1,at=points)
  axis(2,at=c(0,5,10,15))
  
  
  
  
  
  if (generateSVG) {
    dev.off()
  }
  
  
  
}
