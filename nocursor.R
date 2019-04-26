
source('shared.R')

# FIGURE -----

plotReachAftereffects <- function(generateSVG=FALSE, selectPerformance=TRUE, addNormalFits=FALSE) {
  
  if (generateSVG) {
    installed.list <- rownames(installed.packages())
    if ('svglite' %in% installed.list) {
      library('svglite')
      svglite(file='Fig4.svg', width=7.5, height=3, system_fonts=list(sans='Arial', mono='Times New Roman'))
    } else {
      generateSVG=FALSE
    }
  }
  
  par(mfrow=c(1,3))
  
  points <- c(15,25,35,45,55,65,75)
  
  exposure_ini <- getReachAftereffects('exposure',part='initial', selectPerformance=selectPerformance)
  exposure_rem <- getReachAftereffects('exposure',part='remainder', selectPerformance=selectPerformance)
  classic  <- getReachAftereffects('classic',part='all', selectPerformance=selectPerformance)
  exposure <- getReachAftereffects('exposure',part='all', selectPerformance=selectPerformance)
  
  exposureAVGini <- aggregate(endpoint_angle ~ target, data=exposure_ini, FUN=mean) # error?
  exposureAVGrem <- aggregate(endpoint_angle ~ target, data=exposure_rem, FUN=mean)
  classicAVG <- aggregate(endpoint_angle ~ target, data=classic, FUN=mean)
  exposureAVG <- aggregate(endpoint_angle ~ target, data=exposure, FUN=mean)

  exposureCIini <- matrix(unlist(by(exposure_ini$endpoint_angle, INDICES=c(exposure_ini$target), FUN=t.interval)),nrow=2)
  exposureCIrem <- matrix(unlist(by(exposure_rem$endpoint_angle, INDICES=c(exposure_rem$target), FUN=t.interval)),nrow=2)
  classicCI <- matrix(unlist(by(classic$endpoint_angle, INDICES=c(classic$target), FUN=t.interval)),nrow=2)
  exposureCI <- matrix(unlist(by(exposure$endpoint_angle, INDICES=c(exposure$target), FUN=t.interval)),nrow=2)
  
  X <- c(points, rev(points))
  expYini <- c(exposureCIini[1,],rev(exposureCIini[2,]))
  expYrem <- c(exposureCIrem[1,],rev(exposureCIrem[2,]))
  claY <- c(classicCI[1,],rev(classicCI[2,]))
  expY <- c(exposureCI[1,],rev(exposureCI[2,]))
  
  plot(-1000,-1000, main='reach aftereffect decay', xlab='target angle [°]', ylab='reach endpoint deviation [°]', xlim=c(10,80), ylim=c(0,15), axes=F)
  
  mtext('A', side=3, outer=TRUE, at=c(0,1), line=-1, adj=0, padj=1)
  
  polygon(X,expYrem,border=NA,col=colorset[['expActT']])
  polygon(X,expYini,border=NA,col=colorset[['expPasT']])
  
  lines(points,exposureAVGrem$endpoint_angle,col=colorset[['expActS']],lty=2,lwd=1.5)
  lines(points,exposureAVGini$endpoint_angle,col=colorset[['expPasS']],lwd=1.5)
  
  axis(1,at=points)
  axis(2,at=c(0,5,10,15))
  
  legend(10,15,c('iteration 1','iterations 2-5'),col=c(colorset[['expPasS']],colorset[['expActS']]),lty=c(1,2),lwd=c(1.5,1.5),bty='n')
  
  plot(-1000,-1000, main='reach aftereffects', xlab='target angle [°]', ylab='reach endpoint deviation [°]', xlim=c(10,80), ylim=c(0,15), axes=F)
  
  mtext('B', side=3, outer=TRUE, at=c(1/3,1), line=-1, adj=0, padj=1)
  
  polygon(X,claY, border=NA,col=colorset[['claActT']])
  polygon(X,expY, border=NA,col=colorset[['expActT']])
  
  lines(points,classicAVG$endpoint_angle,col=colorset[['claActS']],lwd=1.5)
  lines(points,exposureAVG$endpoint_angle,col=colorset[['expActS']],lwd=1.5)

  axis(1,at=points)
  axis(2,at=c(0,5,10,15))
  
  legend(10,15,c('exposure','classic'),col=c(colorset[['expActS']],colorset[['claActS']]),lty=c(1,1),lwd=c(1.5,1.5),bty='n')
  
  if (addNormalFits) {
    
    points=c(15,25,35,45,55,65,75)
    
    classic  <- getReachAftereffects('classic',part='all', selectPerformance=selectPerformance)
    exposure <- getReachAftereffects('exposure',part='all', selectPerformance=selectPerformance)
    
    # fitting group data:
    exp.fit <- getGaussianFit(x=exposure$target,exposure$endpoint_angle,mu=50,sigma=30,scale=50,offset=4)
    cla.fit <- getGaussianFit(x=classic$target,classic$endpoint_angle,mu=50,sigma=30,scale=50,offset=4)
    
    # get confidence intervals for the peak of the generalization curve for localization shifts:
    cla.RAEshift <- getPeakConfidenceInterval('classic', part='all', CIs=c(.95), selectPerformance=selectPerformance)
    exp.RAEshift <- getPeakConfidenceInterval('exposure', part='all', CIs=c(.95), selectPerformance=selectPerformance)
    
    plot(-1000,-1000, main='generalization curves', xlab='target angle [°]', ylab='reach endpoint deviation [°]', xlim=c(10,80), ylim=c(0,15), axes=F)
    
    mtext('C', side=3, outer=TRUE, at=c(2/3,1), line=-1, adj=0, padj=1)
    
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
    
  }
  
  if (generateSVG) {
    dev.off()
  }
  
}

# VALIDATION (NOT USED ANYMORE) -----

validateReachAftereffects <- function(group, method='45') {
  
  # method can be: '45', 'all_aov', 'all_ttest', 'all_any5deg'
  
  raw.df <- load.DownloadDataframe(url=nocursorURLs[group],filename=sprintf('%s_nocursor.csv',group))
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

exposureNoCursorChange <- function(printChange=FALSE,LMEmethod='chi-squared',selectPerformance=selectPerformance) {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  nc.exp <- getReachAftereffects('exposure', part='initial', difference=FALSE, selectPerformance=selectPerformance)
  
  nc.exp$participant <- factor(nc.exp$participant)
  nc.exp$rotated <- factor(nc.exp$rotated)
  nc.exp$target <- factor(nc.exp$target)
  
  attach(nc.exp)
  
  cat('\nLME with session and target as fixed effects, and participant as random effect:\n\n')
  if (LMEmethod=='chi-squared') {
    print(Anova(lme(endpoint_angle ~ rotated * target, random = ~1|participant, na.action=na.exclude), type=3))
  }
  if (LMEmethod=='Satterthwaite') {
    exp_model_lmer <- lmer(endpoint_angle ~ rotated * target - (1|participant), na.action=na.exclude)
    print(anova(exp_model_lmer,ddf='Satterthwaite',type=3))
  }
  
  detach(nc.exp)

  options('contrasts' <- default.contrasts)
  
  if (printChange) {
    
    nc.exp <- getReachAftereffects('exposure', part='initial', difference=TRUE, selectPerformance=selectPerformance)
    
    nc.exp <- aggregate(endpoint_angle ~ target, data=nc.exp, FUN=mean)
    
    cat(sprintf('\nEXPOSURE\nmax reach aftereffects: %0.1f (at: %d deg)\n',max(nc.exp$endpoint_angle), nc.exp$target[which(nc.exp$endpoint_angle == max(nc.exp$endpoint_angle))]))
    cat(sprintf('avg reach aftereffects: %0.1f (all targets)\n',mean(nc.exp$endpoint_angle, na.rm=TRUE)))
    cat(sprintf('avg reach aftereffects: %0.1f (targets above 15 degrees)\n',mean(nc.exp$endpoint_angle[which(nc.exp$target > 15)], na.rm=TRUE)))
    
    nc.cla <- getReachAftereffects('classic', part='initial', difference=TRUE, selectPerformance=selectPerformance)
    
    nc.cla <- aggregate(endpoint_angle ~ target, data=nc.cla, FUN=mean)
    
    cat(sprintf('\nCLASSIC\nmax reach aftereffects: %0.1f (at: %d deg)\n',max(nc.cla$endpoint_angle), nc.cla$target[which(nc.cla$endpoint_angle == max(nc.cla$endpoint_angle))]))
    cat(sprintf('avg reach aftereffects: %0.1f (all targets)\n',mean(nc.cla$endpoint_angle, na.rm=TRUE)))
    cat(sprintf('avg reach aftereffects: %0.1f (targets above 15 degrees)\n',mean(nc.cla$endpoint_angle[which(nc.cla$target > 15)], na.rm=TRUE)))
    
  }
  
}

exposureAftereffectsPersistent <- function(LMEmethod='chi-squared', selectPerformance=selectPerformance) {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  exp1 <- getReachAftereffects('exposure', part='initial', difference=TRUE, selectPerformance=selectPerformance)
  exp2 <- getReachAftereffects('exposure', part='remainder', difference=TRUE, selectPerformance=selectPerformance)
  
  exp1$iteration <- 1
  exp2$iteration <- 2
  
  exp <- rbind(exp1, exp2)
  
  exp$participant <- factor(exp$participant)
  exp$iteration <- factor(exp$iteration)
  exp$target <- factor(exp$target)
  
  attach(exp)
  
  cat('\nLME with session and target as fixed effects, and participant as random effect:\n\n')
  
  
  if (LMEmethod=='chi-squared') {
    print(Anova(lme(endpoint_angle ~ iteration * target, random = ~1|participant, na.action=na.exclude), type=3))
  }
  if (LMEmethod=='Satterthwaite') {
    exp_model_lmer <- lmer(endpoint_angle ~ iteration * target - (1|participant), na.action=na.exclude)
    print(anova(exp_model_lmer,ddf='Satterthwaite',type=3))
  }
  
  
  detach(exp)
  
  options('contrasts' <- default.contrasts)
  
}

exposureClassicReachAftereffects <- function(noTarget=FALSE, remove15=FALSE, LMEmethod='chi-squared', selectPerformance=selectPerformance) {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  exp <- getReachAftereffects('exposure', part='all', difference=TRUE, selectPerformance=selectPerformance)
  cla <- getReachAftereffects('classic', part='all', difference=TRUE, selectPerformance=selectPerformance)
  if (remove15) {
    exp <- exp[-which(exp$target == 15),]
    cla <- cla[-which(cla$target == 15),]
  }
  
  exp$training <- 1
  cla$training <- 2
  
  RAE <- rbind(exp, cla)
  
  RAE$participant <- factor(RAE$participant)
  RAE$training <- factor(RAE$training)
  RAE$target <- factor(RAE$target)
  
  attach(RAE)
  
  if (noTarget) {
    
    cat('\nLME with training type as fixed effects (without interacting), and participant and target as random effect:\n\n')
    if (LMEmethod=='chi-squared') {
      print(Anova(lme(endpoint_angle ~ training, random = ~1|participant/target, na.action=na.exclude), type=3))
    }
    if (LMEmethod=='Satterthwiate') {
      exp_cla_model_lmer <- lmer(endpoint_angle ~ training - (1|participant/target), na.action=na.exclude)
      print(anova(exp_cla_model_lmer,ddf='Satterthwaite',type=3))
    }
    
  } else {
    
    cat('\nLME with training type and target as fixed effects, and participant as random effect:\n\n')
    
    if (LMEmethod=='chi-squared') {
      print(Anova(lme(endpoint_angle ~ training * target, random = ~1|participant, na.action=na.exclude), type=3))
    }
    if (LMEmethod=='Satterthwaite') {
      exp_cla_model_lmer <- lmer(endpoint_angle ~ training * target - (1|participant), na.action=na.exclude)
      print(anova(exp_cla_model_lmer,ddf='Satterthwaite',type=3))
    }
    
  }
  
  detach(RAE)
  
  options('contrasts' <- default.contrasts)
  
}

# BLINK DETECTION PERFORMANCE -----

plotBlinkDetection <- function(generateSVG=FALSE) {
  
  if (generateSVG) {
    installed.list <- rownames(installed.packages())
    if ('svglite' %in% installed.list) {
      library('svglite')
      svglite(file='Fig2.5.svg', width=4, height=4, system_fonts=list(sans='Arial', mono='Times New Roman'))
    } else {
      generateSVG=FALSE
    }
  }
  
  #blinks <- read.csv(, stringsAsFactors=F)
  blinks <- load.DownloadDataframe(informationURLs['blinkdetect'],'blinkdetect_exposure.csv')
  
  participants <- unique(blinks$participant)
  plot(x=c(1,2),y=c(.5,.5),type='l',lty=2,col=rgb(0,0,0),xlim=c(0,3),ylim=c(0,1),main='cursor-blink detection',axes=F,xlab='task',ylab='performance')
  
  for(participant in participants) {
    pblink <- blinks[which(blinks$participant == participant),]
    lines(x=c(1:2),y=pblink$performance,col=rgb(0.7,0.75,0.7))
  }
  
  blinks$taskno <- blinks$taskno + (blinks$rotated_b)
  avg <- aggregate(performance ~ taskno, data=blinks, FUN=mean)
  lines(x=c(1:2),y=avg$performance,col=rgb(0,0,0))
  
  axis(side=1, at=c(1:2), labels=c('aligned','rotated'))
  axis(side=2, at=c(0,0.5,1))
  
  if (generateSVG) {
    dev.off()
  }
  
}

# GENERALIZATION PEAK CI -----

getPeakConfidenceInterval <- function(group,validate=FALSE,part='all',CIs=c(.95), selectPerformance=selectPerformance,iterations=1000) {
  
  filename <- sprintf('maxima_RAE_%s.csv', group)
  
  if (file.exists(filename)) {
    
    #cat(sprintf('\nloading peak RAE generalization from file for: %s\n',toupper(group)))
    
    df <- read.csv(filename, stringsAsFactors=FALSE)
    
  } else {
    
    cat(sprintf('\nbootstrapping peak RAE generalization for: %s\n',toupper(group)))
    
    RAE <- getReachAftereffects(group, part=part, selectPerformance=selectPerformance)
    
    # if (validate) {
    #   validParticipants <- validateReachAftereffects(group, method='45')
    #   validParticipants <- validParticipants$participant[which(validParticipants$RAE == 1)]
    #   RAE <- RAE[which(RAE$participant %in% validParticipants),]
    # }
    
    RAE <- xtabs(endpoint_angle~.,RAE)
    
    # data <- bootstrapGaussianPeak(data=RAE,bootstraps=iterations,mu=47.5,sigma=30,scale=10,offset=4,CIs=CIs)
    data <- bootstrapGaussianPeak(data=RAE,bootstraps=iterations,mu=50,sigma=30,scale=50,offset=4,CIs=CIs)
    
    df <- data.frame('level'=names(data),'value'=data)
    
    write.csv(df,filename,row.names=FALSE,quote=FALSE)
    
  }
  
  return(df)
  
}

# DATA DESCRIPTIVES -----

countSelectedNoCursors <- function(group, ignoreRepetitions=FALSE, selectPerformance=TRUE) {
  
  df <- load.DownloadDataframe(url=nocursorURLs[group],filename=sprintf('nocursor_%s.csv',group))
  
  if (selectPerformance & group=='exposure') {
    blinks <- load.DownloadDataframe(informationURLs['blinkdetect'],'blinkdetect_exposure.csv')
    OKparticipants <- blinks$participant[which(blinks$rotated_b == 1 & blinks$performance > 0.65)]
    df <- df[which(df$participant %in% OKparticipants),]
  }
  
  participant <- c()
  rotated <- c()
  passive <- c()
  repetition <- c()
  trials <- c()
  
  mintrials <- 21
  
  participants <- unique(df$participant)
  
  for (ppid in participants) {
    
    ppdf <- df[which(df$participant == ppid),]
    
    for (session in c(0,1)) {
      
      subdf <- ppdf[which(ppdf$rotated == session),]
      
      iters <- unique(subdf$repetition)
      
      for (iterno in c(1:length(iters))) {
        
        iter <- iters[iterno]
        
        iterdf <- subdf[which(subdf$repetition == iter),]
        
        Ntrials <- dim(iterdf)[1]
        
        if (Ntrials < mintrials) {
          mintrials <- Ntrials
        }
        
        participant <- c(participant, ppid)
        rotated     <- c(rotated, session)
        repetition  <- c(repetition, iter)
        trials      <- c(trials, (Ntrials/.21))
        
      }
      
    }
    
  }
  
  #cat(sprintf('\nminimum trials selected: %d\n\n',mintrials))
  
  return(data.frame(participant, rotated, repetition, trials))
  
}