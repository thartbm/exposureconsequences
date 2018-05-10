
source('shared.R')

# first make sure we have all necessary packages installed and loaded:
# installRequire.Packages(c())

# not yet sure if we'll use Chi-square, or Satterthwaite estimates of F/p-values:
# installRequire.Packages(c('nlme', 'car', 'lme4', 'lmerTest'))

required.packages = c('nlme', 'car', 'lme4', 'lmerTest', 'svglite')
installRequire.Packages(required.packages)

# nocursor data can be downloaded from OSF:
# groupURLs <- c('exposure'='https://osf.io/9s6au/?action=download', 'classic'='https://osf.io/8hm7f/?action=download')


# # plot both groups reach aftereffects in one figure
# OLDplotReachAftereffects <- function(validate=FALSE) {
#   
#   points <- c(15,25,35,45,55,65,75)
#   
#   exposure <- getReachAftereffects('exposure',part='initial')
#   classic  <- getReachAftereffects('classic',part='initial')
#   exposure_rem <- getReachAftereffects('exposure',part='remainder')
#   classic_rem  <- getReachAftereffects('classic',part='remainder')
#   
#   if (validate) {
#     
#     validExp <- validateReachAftereffects('exposure')
#     validExpParticipants <- validExp$participant[which(validExp$RAE == 1)]
#     exposure <- exposure[which(exposure$participant %in% validExpParticipants),]
# 
#     validCla <- validateReachAftereffects('classic')
#     validClaParticipants <- validCla$participant[which(validCla$RAE == 1)]
#     classic <- classic[which(classic$participant %in% validClaParticipants),]
#     
#   }
#   
#   exposureAVG <- aggregate(endpoint_angle ~ target, data=exposure, FUN=mean)
#   classicAVG <- aggregate(endpoint_angle ~ target, data=classic, FUN=mean)
#   
#   exposureAVGrem <- aggregate(endpoint_angle ~ target, data=exposure_rem, FUN=mean)
#   classicAVGrem <- aggregate(endpoint_angle ~ target, data=classic_rem, FUN=mean)
#   
#   exposureCI <- matrix(unlist(by(exposure$endpoint_angle, INDICES=c(exposure$target), FUN=t.interval)),nrow=2)
#   classicCI <- matrix(unlist(by(classic$endpoint_angle, INDICES=c(classic$target), FUN=t.interval)),nrow=2)
#   
#   X <- c(points, rev(points))
#   expY <- c(exposureCI[1,],rev(exposureCI[2,]))
#   claY <- c(classicCI[1,],rev(classicCI[2,]))
#   
#   plot(-1000,-1000, main='reach aftereffects', xlab='target angle [deg]', ylab='reach endpoint deviation [deg]', xlim=c(10,80), ylim=c(0,15), axes=F)
#   
#   polygon(X,claY,border=NA,col=rgb(1,0,0,.2))
#   polygon(X,expY,border=NA,col=rgb(0,0,1,.2))
# 
#   lines(points,classicAVG$endpoint_angle,col=rgb(1,0,0))
#   lines(points,exposureAVG$endpoint_angle,col=rgb(0,0,1))
#   lines(points,classicAVGrem$endpoint_angle,col=rgb(1,0,0),lty=2)
#   lines(points,exposureAVGrem$endpoint_angle,col=rgb(0,0,1),lty=2)
# 
#   axis(1,at=points)
#   axis(2,at=c(0,5,10,15))
#   
#   legend(10,15,c('exposure','classic'),col=c(rgb(0,0,1),rgb(1,0,0)),lty=c(1,1),bty='n')
#   
# }

plotReachAftereffects <- function(generateSVG=FALSE) {
  
  if (generateSVG) {
    installed.list <- rownames(installed.packages())
    if ('svglite' %in% installed.list) {
      svglite(file='Fig1.svg', width=7.5, height=2.5, system_fonts=list(sans='Arial', mono='Times New Roman'))
    } else {
      generateSVG=FALSE
    }
  }
  
  par(mfrow=c(1,3))
  
  points <- c(15,25,35,45,55,65,75)
  
  exposure_ini <- getReachAftereffects('exposure',part='initial')
  exposure_rem <- getReachAftereffects('exposure',part='remainder')
  classic  <- getReachAftereffects('classic',part='all')
  exposure <- getReachAftereffects('exposure',part='all')
  
  exposureAVGini <- aggregate(endpoint_angle ~ target, data=exposure_ini, FUN=mean)
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
  
  plot(-1000,-1000, main='decay of reach aftereffects', xlab='target angle [deg]', ylab='reach endpoint deviation [deg]', xlim=c(10,80), ylim=c(0,15), axes=F)
  
  polygon(X,expYrem,border=NA,col=colorset[['expActT']])
  polygon(X,expYini,border=NA,col=colorset[['expPasT']])
  
  lines(points,exposureAVGrem$endpoint_angle,col=colorset[['expActS']],lty=2,lwd=2)
  lines(points,exposureAVGini$endpoint_angle,col=colorset[['expPasS']],lwd=2)
  
  axis(1,at=points)
  axis(2,at=c(0,5,10,15))
  
  legend(10,15,c('immediate (iteration 1)','delayed (iterations 2-5)'),col=c(colorset[['expPasS']],colorset[['expActS']]),lty=c(1,2),lwd=c(2,2),bty='n')
  
  plot(-1000,-1000, main='reach aftereffects', xlab='target angle [deg]', ylab='reach endpoint deviation [deg]', xlim=c(10,80), ylim=c(0,15), axes=F)
  
  polygon(X,claY, border=NA,col=colorset[['claActT']])
  polygon(X,expY, border=NA,col=colorset[['expActT']])
  
  lines(points,classicAVG$endpoint_angle,col=colorset[['claActS']],lwd=2)
  lines(points,exposureAVG$endpoint_angle,col=colorset[['expActS']],lwd=2)

  axis(1,at=points)
  axis(2,at=c(0,5,10,15))
  
  legend(10,15,c('exposure','classic'),col=c(colorset[['expActS']],colorset[['claActS']]),lty=c(1,1),lwd=c(2,2),bty='n')
  
  if (generateSVG) {
    dev.off()
  }
  
}

getPeakConfidenceInterval <- function(group,validate=FALSE,part='initial',CIs=c(.95)) {
  
  cat(sprintf('\n%s\n\n',toupper(group)))
  
  RAE <- getReachAftereffects(group, part=part)
  
  if (validate) {
    validParticipants <- validateReachAftereffects(group, method='45')
    validParticipants <- validParticipants$participant[which(validParticipants$RAE == 1)]
    RAE <- RAE[which(RAE$participant %in% validParticipants),]
  }
  
  RAE <- xtabs(endpoint_angle~.,RAE)
  
  bootstrapGaussianPeak(data=RAE,bootstraps=1000,mu=47.5,sigma=30,scale=10,offset=4,CIs=CIs)
  
}

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

exposureNoCursorChange <- function(printChange=FALSE,LMEmethod='chi-squared') {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  nc.exp <- getReachAftereffects('exposure', part='initial', difference=FALSE)
  
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
    
    nc.exp <- getReachAftereffects('exposure', part='initial', difference=TRUE)
    
    nc.exp <- aggregate(endpoint_angle ~ target, data=nc.exp, FUN=mean)
    
    cat(sprintf('\nEXPOSURE\nmax reach aftereffects: %0.1f (at: %d deg)\n',max(nc.exp$endpoint_angle), nc.exp$target[which(nc.exp$endpoint_angle == max(nc.exp$endpoint_angle))]))
    cat(sprintf('avg reach aftereffects: %0.1f (all targets)\n',mean(nc.exp$endpoint_angle, na.rm=TRUE)))
    cat(sprintf('avg reach aftereffects: %0.1f (targets above 15 degrees)\n',mean(nc.exp$endpoint_angle[which(nc.exp$target > 15)], na.rm=TRUE)))
    
    nc.cla <- getReachAftereffects('classic', part='initial', difference=TRUE)
    
    nc.cla <- aggregate(endpoint_angle ~ target, data=nc.cla, FUN=mean)
    
    cat(sprintf('\nCLASSIC\nmax reach aftereffects: %0.1f (at: %d deg)\n',max(nc.cla$endpoint_angle), nc.cla$target[which(nc.cla$endpoint_angle == max(nc.cla$endpoint_angle))]))
    cat(sprintf('avg reach aftereffects: %0.1f (all targets)\n',mean(nc.cla$endpoint_angle, na.rm=TRUE)))
    cat(sprintf('avg reach aftereffects: %0.1f (targets above 15 degrees)\n',mean(nc.cla$endpoint_angle[which(nc.cla$target > 15)], na.rm=TRUE)))
    
  }
  
}

exposureAftereffectsPersistent <- function(LMEmethod='chi-squared') {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  exp1 <- getReachAftereffects('exposure', part='initial', difference=TRUE)
  exp2 <- getReachAftereffects('exposure', part='remainder', difference=TRUE)
  
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

exposureClassicReachAftereffects <- function(noTarget=FALSE, remove15=FALSE, LMEmethod='chi-squared') {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  exp <- getReachAftereffects('exposure', part='all', difference=TRUE)
  cla <- getReachAftereffects('classic', part='all', difference=TRUE)
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
      exp_cla_model_lmer <- lmer(endpoint_angle ~ training - (1|participant), na.action=na.exclude)
      print(anova(exp_cla_model_lmer,ddf='Satterthwaite',type=3))
    }
    
  }
  
  detach(RAE)
  
  options('contrasts' <- default.contrasts)
  
}
