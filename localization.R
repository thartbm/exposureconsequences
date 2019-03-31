
source('shared.R')


# PLOT / FIGURE ------

plotLocalization <- function(classicOnline=FALSE, generateSVG=FALSE, selectPerformance=TRUE, remove15=FALSE, thirdPanel='2x2') {
  
  # thirdPanel=
  # 'peakCIs'
  # 'classicOnline'
  # '2x2'
  
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
  
  # get the confidence intervals for polygon areas:
  exp.act <- exp[which(exp$passive_b == 0),]
  exp.CI.act <- matrix(unlist(by(exp.act$taperror_deg, INDICES=c(exp.act$handangle_deg), FUN=t.interval)),nrow=2)
  exp.pas <- exp[which(exp$passive_b == 1),]
  exp.CI.pas <- matrix(unlist(by(exp.pas$taperror_deg, INDICES=c(exp.pas$handangle_deg), FUN=t.interval)),nrow=2)
  cla.act <- cla[which(cla$passive_b == 0),]
  cla.CI.act <- matrix(unlist(by(cla.act$taperror_deg, INDICES=c(cla.act$handangle_deg), FUN=t.interval)),nrow=2)
  cla.pas <- cla[which(cla$passive_b == 1),]
  cla.CI.pas <- matrix(unlist(by(cla.pas$taperror_deg, INDICES=c(cla.pas$handangle_deg), FUN=t.interval)),nrow=2)

  if (thirdPanel == 'classicOnline') {
    onl <- getPointLocalization('online', difference=TRUE, verbose=FALSE, selectPerformance=FALSE)
    
    onl.avg <- aggregate(taperror_deg ~ passive_b + handangle_deg, data=onl, FUN=mean)
    onl.avg.act <- onl.avg[which(onl.avg$passive_b == 0),]
    onl.avg.pas <- onl.avg[which(onl.avg$passive_b == 1),]
    
    onl.act <- onl[which(onl$passive_b == 0),]
    onl.CI.act <- matrix(unlist(by(onl.act$taperror_deg, INDICES=c(onl.act$handangle_deg), FUN=t.interval)),nrow=2)
    onl.pas <- onl[which(onl$passive_b == 1),]
    onl.CI.pas <- matrix(unlist(by(onl.pas$taperror_deg, INDICES=c(onl.pas$handangle_deg), FUN=t.interval)),nrow=2)
  }
  
  if (generateSVG) {
    installed.list <- rownames(installed.packages())
    if ('svglite' %in% installed.list) {
      library('svglite')
      svglite(file='Fig3.svg', width=7.5, height=3, system_fonts=list(sans='Arial', mono='Times New Roman'))
    } else {
      generateSVG=FALSE
    }
  }
  
  par(mfrow=c(1,3))
  
  points <- c(15,25,35,45,55,65,75)
  
  # panel A: exposure localization (active vs. passive)
  
  plot(-1000,-1000, main='exposure', xlab='hand angle [°]', ylab='localization shift [°]', xlim=c(10,80), ylim=c(0,-15), axes=F)
  
  mtext('A', side=3, outer=TRUE, at=c(0,1), line=-1, adj=0, padj=1)
  
  X <- c(points, rev(points))
  exp.act.Y <- c(exp.CI.act[1,],rev(exp.CI.act[2,]))
  exp.pas.Y <- c(exp.CI.pas[1,],rev(exp.CI.pas[2,]))
  
  polygon(X,exp.act.Y,border=NA,col=colorset[['expActT']])
  polygon(X,exp.pas.Y,border=NA,col=colorset[['expPasT']])
  
  lines(points[1:2],exp.avg.act$taperror_deg[1:2],col=colorset[['expActS']],lty=2,lwd=1.5)
  lines(points[2:7],exp.avg.act$taperror_deg[2:7],col=colorset[['expActS']],lty=1,lwd=1.5)
  lines(points[1:2],exp.avg.pas$taperror_deg[1:2],col=colorset[['expPasS']],lty=2,lwd=1.5)
  lines(points[2:7],exp.avg.pas$taperror_deg[2:7],col=colorset[['expPasS']],lty=1,lwd=1.5)
  
  axis(1,at=points)
  axis(2,at=c(0,-5,-10,-15))
  
  legend(10,-15,c('passive','active'),col=c(colorset[['expPasS']],colorset[['expActS']]),lty=c(1,1),lwd=c(1.5,1.5),bty='n')
  
  # panel B: classic localization (active vs. passive)
  
  plot(-1000,-1000, main='classic', xlab='hand angle [°]', ylab='localization shift [°]', xlim=c(10,80), ylim=c(0,-15), axes=F)
  
  mtext('B', side=3, outer=TRUE, at=c(1/3,1), line=-1, adj=0, padj=1)
  
  X <- c(points, rev(points))
  cla.act.Y <- c(cla.CI.act[1,],rev(cla.CI.act[2,]))
  cla.pas.Y <- c(cla.CI.pas[1,],rev(cla.CI.pas[2,]))
  
  polygon(X,cla.act.Y,border=NA,col=colorset[['claActT']])
  polygon(X,cla.pas.Y,border=NA,col=colorset[['claPasT']])
  
  lines(points[1:2],cla.avg.act$taperror_deg[1:2],col=colorset[['claActS']],lty=2,lwd=1.5)
  lines(points[2:7],cla.avg.act$taperror_deg[2:7],col=colorset[['claActS']],lty=1,lwd=1.5)
  lines(points[1:2],cla.avg.pas$taperror_deg[1:2],col=colorset[['claPasS']],lty=2,lwd=1.5)
  lines(points[2:7],cla.avg.pas$taperror_deg[2:7],col=colorset[['claPasS']],lty=1,lwd=1.5)
  
  axis(1,at=points)
  axis(2,at=c(0,-5,-10,-15))
  
  legend(10,-15,c('passive','active'),col=c(colorset[['claPasS']],colorset[['claActS']]),lty=c(1,1),lwd=c(1.5,1.5),bty='n')
  
  if (thirdPanel == 'classicOnline') {
    
    plot(-1000,-1000, main='online', xlab='hand angle [°]', ylab='localization shift [°]', xlim=c(10,80), ylim=c(0,-15), axes=F)
    
    mtext('C', side=3, outer=TRUE, at=c(2/3,1), line=-1, adj=0, padj=1)
    
    X <- c(points, rev(points))
    onl.act.Y <- c(onl.CI.act[1,],rev(onl.CI.act[2,]))
    onl.pas.Y <- c(onl.CI.pas[1,],rev(onl.CI.pas[2,]))
    
    polygon(X,onl.act.Y,border=NA,col=colorset[['onlActT']])
    polygon(X,onl.pas.Y,border=NA,col=colorset[['onlPasT']])
    
    lines(points[1:2],onl.avg.act$taperror_deg[1:2],col=colorset[['onlActS']],lty=2,lwd=1.5)
    lines(points[2:7],onl.avg.act$taperror_deg[2:7],col=colorset[['onlActS']],lty=1,lwd=1.5)
    lines(points[1:2],onl.avg.pas$taperror_deg[1:2],col=colorset[['onlPasS']],lty=2,lwd=1.5)
    lines(points[2:7],onl.avg.pas$taperror_deg[2:7],col=colorset[['onlPasS']],lty=1,lwd=1.5)
    
    axis(1,at=points)
    axis(2,at=c(0,-5,-10,-15))
    
    legend(10,-15,c('passive','active'),col=c(colorset[['onlPasS']],colorset[['onlActS']]),lty=c(1,1),lwd=c(1.5,1.5),bty='n')
    
  } else if (thirdPanel == 'peakCIs') {
    
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
    
    plot(-1000,-1000, main='generalization curves', xlab='hand angle [°]', ylab='localization shift [°]', xlim=c(10,80), ylim=c(0,-15), axes=F)
    
    mtext('C', side=3, outer=TRUE, at=c(2/3,1), line=-1, adj=0, padj=1)
    
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
    
    # add tick marks:
    axis(1,at=points)
    axis(2,at=c(0,-5,-10,-15))
    
  } else if (thirdPanel == '2x2') {
    
    # get the data again, because we remove the 15 degree point?
    #exp <- getPointLocalization('exposure', difference=TRUE, verbose=FALSE, selectPerformance=selectPerformance)
    #cla <- getPointLocalization('classic', difference=TRUE, verbose=FALSE, selectPerformance=FALSE)
    
    # get the averages for the line plots:
    exp.act <- aggregate(taperror_deg ~ participant, data=exp[which(exp$passive_b == 0 & exp$handangle_deg > 15),], FUN=mean)
    exp.pas <- aggregate(taperror_deg ~ participant, data=exp[which(exp$passive_b == 1 & exp$handangle_deg > 15),], FUN=mean)
    #
    cla.act <- aggregate(taperror_deg ~ participant, data=cla[which(cla$passive_b == 0 & exp$handangle_deg > 15),], FUN=mean)
    cla.pas <- aggregate(taperror_deg ~ participant, data=cla[which(cla$passive_b == 1 & exp$handangle_deg > 15),], FUN=mean)
    
    # get the confidence intervals for polygon areas:
    exp.act.CI <- t.interval(exp.act$taperror_deg)
    exp.pas.CI <- t.interval(exp.pas$taperror_deg)
    cla.act.CI <- t.interval(cla.act$taperror_deg)
    cla.pas.CI <- t.interval(cla.pas$taperror_deg)
    
    plot(-1000,-1000, main='localization shifts', xlab='localization task', ylab='mean localization shift [°]', xlim=c(10,80), ylim=c(0,-15), axes=F)
    
    mtext('C', side=3, outer=TRUE, at=c(2/3,1), line=-1, adj=0, padj=1)
    
    X    <- c(20,70,70,20)
    Yexp <- c(exp.act.CI[1],exp.pas.CI[1],exp.pas.CI[2],exp.act.CI[2]) 
    Ycla <- c(cla.act.CI[1],cla.pas.CI[1],cla.pas.CI[2],cla.act.CI[2]) 
    
    polygon(X,Yexp,border=NA,col=colorset[['expActT']])
    polygon(X,Ycla,border=NA,col=colorset[['claActT']])
    
    lines(c(20,70),c(mean(exp.act$taperror_deg),mean(exp.pas$taperror_deg)),col=colorset[['expActS']],lty=1,lw=1.5)
    lines(c(20,70),c(mean(cla.act$taperror_deg),mean(cla.pas$taperror_deg)),col=colorset[['claActS']],lty=1,lw=1.5)
    
    legend(10,-15,c('exposure','classic'),col=c(colorset[['expActS']],colorset[['claActS']]),lty=c(1,1),lwd=c(1.5,1.5),bty='n')
    
    # add tick marks:
    axis(1,at=c(20,70),labels=c('active','passive'))
    axis(2,at=c(0,-5,-10,-15))
    
  }
  
  if (generateSVG) {
    dev.off()
  }
  
}

plotALignedRotatedLocalization <- function(classicOnline=FALSE, generateSVG=FALSE, selectPerformance=TRUE, remove15=FALSE) {
  
  groups <- c('exposure','classic')
  
  if (classicOnline) {
    
    groups <- c(groups,'online')
    
  }
  
  if (generateSVG) {
    installed.list <- rownames(installed.packages())
    if ('svglite' %in% installed.list) {
      # nothing
      library('svglite')
    } else {
      generateSVG=FALSE
    }
  }
  
  par(mfrow=c(7,3),mai=c(.6,.5,.01,.01))
  
  points <- c(15,25,35,45,55,65,75)
  
  
  for (group in groups) {
    
    SP <- FALSE
    if (group == 'exposure') {
      SP <- selectPerformance
    }
    
    df <- load.DownloadDataframe(url=localizationURLs[group],filename=sprintf('localization_%s.csv',group))
    
    if (selectPerformance & group=='exposure') {
      blinks <- load.DownloadDataframe(informationURLs['blinkdetect'],'blinkdetect_exposure.csv')
      OKparticipants <- blinks$participant[which(blinks$rotated_b == 1 & blinks$performance > 0.65)]
      df <- df[which(df$participant %in% OKparticipants),]
    }
    
    participants <- unique(df$participant)
    
    for (passive in c(0,1)) {
      
      if (generateSVG) {
        # should have been loaded if available:
        svglite(file=sprintf('Fig7_%s_%s.svg',group,c('active','passive')[passive+1]), width=8.5, height=11, system_fonts=list(sans='Arial', mono='Times New Roman'))
        par(mfrow=c(7,3),mai=c(.6,.5,.01,.01))
      }
        
      for (participant in participants) {
      
        # create plot:
        plot(-1000,-1000, main='', xlab='hand angle [°]', ylab='localization error [°]', xlim=c(0,90), ylim=c(20,-40), axes=F)
        
        lines(c(10,80),c(0,0),col='#999999')
        
        for (rotated in c(0,1)) {
        
          subdf <- df[which(df$participant == participant & df$rotated_b == rotated & df$passive_b == passive),]
          
          color <- colorset[[sprintf('%s%s%s',substr(group,1,3),c('Act','Pas')[passive+1],c('T','S')[rotated+1])]]
          
          points(subdf$handangle_deg, subdf$taperror_deg,col=color)
          
          if (nrow(subdf) == 0) {
            # participant has no data in sub-condition?
            next()
          }
          
          locdf <- getLocalizationPoints(subdf, points=points, removeOutliers=TRUE)
          
          lty <- c(2,1)[rotated]
          lines(locdf$handangle_deg,locdf$taperror_deg,lty=1,lw=2,col=color)
          
        }
        
        axis(1,at=points)
        axis(2,at=c(10,-10,-30))
        
      }
      
      if (generateSVG) {
        dev.off()
      }
        
    }
    
  }
  
}

# ANALYSES ------

exposureLocalization <- function(remove15=TRUE, LMEmethod='chi-squared', selectPerformance=TRUE) {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  exp <- getPointLocalization('exposure', difference=FALSE, verbose=FALSE, selectPerformance=selectPerformance)
  
  if (remove15) {
    exp <- exp[-which(exp$handangle_deg == 15),]
  }
  
  exp$participant   <- factor(exp$participant)
  exp$rotated_b     <- factor(exp$rotated_b)
  exp$passive_b     <- factor(exp$passive_b)
  exp$handangle_deg <- factor(exp$handangle_deg)
  
  attach(exp)
  
  cat('\nLME with session, target and movement type as fixed effects, and participant as random effect:\n\n')
  
  if (LMEmethod=='chi-squared') {
    print(Anova(lme(taperror_deg ~ rotated_b * passive_b * handangle_deg, random = ~1|participant, na.action=na.exclude), type=3))
  }
  if (LMEmethod=='Satterthwaite') {
    exp_model_lmer <- lmer(taperror_deg ~ rotated_b * passive_b * handangle_deg - (1|participant), na.action=na.exclude)
    print(anova(exp_model_lmer,ddf='Satterthwaite',type=3))
  }
  
  detach(exp)
  
  options('contrasts' <- default.contrasts)
  
}

exposureLocalizationShift <- function(noHandAngle=FALSE, remove15=TRUE, LMEmethod='chi-squared', selectPerformance=TRUE) {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  exp <- getPointLocalization('exposure', difference=TRUE, verbose=FALSE, selectPerformance=selectPerformance)
  
  if (remove15) {
    exp <- exp[-which(exp$handangle_deg == 15),]
  }
  
  exp$participant   <- factor(exp$participant)
  exp$passive_b     <- factor(exp$passive_b)
  exp$handangle_deg <- factor(exp$handangle_deg)
  
  attach(exp)
  
  if (noHandAngle) {
    
    cat('\nLME with movement type as fixed effects - ignoring hand angle, and participant as random effect:\n\n')
    
    if (LMEmethod=='chi-squared') {
      print(Anova(lme(taperror_deg ~ passive_b, random = ~1|participant, na.action=na.exclude), type=3))
    }
    if (LMEmethod=='Satterthwaite') {
      exp_model_lmer <- lmer(taperror_deg ~ passive_b - (1|participant), na.action=na.exclude)
      print(anova(exp_model_lmer,ddf='Satterthwaite',type=3))
    }
    
  } else {
    
    cat('\nLME with hand angle and movement type as fixed effects, and participant as random effect:\n\n')
    
    if (LMEmethod=='chi-squared') {
      print(Anova(lme(taperror_deg ~ passive_b * handangle_deg, random = ~1|participant, na.action=na.exclude), type=3))
    }
    if (LMEmethod=='Satterthwaite') {
      exp_model_lmer <- lmer(taperror_deg ~ passive_b * handangle_deg - (1|participant), na.action=na.exclude)
      print(anova(exp_model_lmer,ddf='Satterthwaite',type=3))
    }
    
  }
  
  detach(exp)
  
  options('contrasts' <- default.contrasts)
  
}

groupLocalization <- function(model='full', remove15=TRUE, LMEmethod='chi-squared', selectPerformance=TRUE) {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  exp <- getPointLocalization('exposure', difference=TRUE, verbose=FALSE, selectPerformance=selectPerformance)
  cla <- getPointLocalization('classic', difference=TRUE, verbose=FALSE, selectPerformance=selectPerformance)
  
  loc <- rbind(exp, cla)
  
  if (remove15) {
    loc <- loc[-which(loc$handangle_deg == 15),]
  }
  
  loc$group         <- factor(loc$group)
  loc$participant   <- factor(loc$participant)
  loc$passive_b     <- factor(loc$passive_b)
  loc$handangle_deg <- factor(loc$handangle_deg)
  
  attach(loc)
  
  if (model == 'full') {
    cat('\nLME with group, hand angle and movement type as fixed effects, and participant as random effect:\n\n')
    if (LMEmethod=='chi-squared') {
      print(Anova(lme(taperror_deg ~ group * passive_b * handangle_deg, random = ~1|participant, na.action=na.exclude), type=3))
    }
    if (LMEmethod=='Satterthwaite') {
      loc_model_lmer <- lmer(taperror_deg ~ group * passive_b * handangle_deg - (1|participant), na.action=na.exclude)
      print(anova(loc_model_lmer,ddf='Satterthwaite',type=3))
    }
  }
  if (model == 'restricted') {
    cat('\nLME with three terms only, removing some main effects and interactions:\n\n')
    if (LMEmethod=='chi-squared') {
      print(Anova(lme(taperror_deg ~ group + group:passive_b + group:handangle_deg, random = c(~1|passive_b, ~1|handangle_deg, ~1|participant), na.action=na.exclude), type=3))
    }
    if (LMEmethod=='Satterthwaite') {
      loc_model_lmer <- lmer(taperror_deg ~ group + group:passive_b + group:handangle_deg - (1|participant), na.action=na.exclude)
      print(anova(loc_model_lmer,ddf='Satterthwaite',type=3))
    }
  }
  if (model == 'handangle') {
    cat('\nLME with group and hand angle as fixed effects, and participant and movement type as random effects:\n\n')
    if (LMEmethod=='chi-squared') {
      print(Anova(lme(taperror_deg ~ group * handangle_deg, random = ~1|participant/taperror_deg, na.action=na.exclude), type=3))
    }
    if (LMEmethod=='Satterthwaite') {
      loc_model_lmer <- lmer(taperror_deg ~ group * handangle_deg - (1|participant), na.action=na.exclude)
      print(anova(loc_model_lmer,ddf='Satterthwaite',type=3))
    }
  }
  if (model == 'movementtype') {
    cat('\nLME with group and movement type as fixed effects and participant and hand angle as random effects:\n\n')
    if (LMEmethod=='chi-squared') {
      print(Anova(lme(taperror_deg ~ group * passive_b, random = ~1|participant/handangle_deg, na.action=na.exclude), type=3))
    }
    if (LMEmethod=='Satterthwaite') {
      loc_model_lmer <- lmer(taperror_deg ~ group * passive_b - (1|participant), na.action=na.exclude)
      print(anova(loc_model_lmer,ddf='Satterthwaite',type=3))
    }
  }
  
  detach(loc)
  
  options('contrasts' <- default.contrasts)
  
}

classicLocalizationShift <- function(factors='both',remove15=TRUE, LMEmethod='chi-squared') {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  cla <- getPointLocalization('classic', difference=TRUE, verbose=FALSE)
  
  if (remove15) {
    cla <- cla[-which(cla$handangle_deg == 15),]
  }
  
  cla$participant   <- factor(cla$participant)
  cla$passive_b     <- factor(cla$passive_b)
  cla$handangle_deg <- factor(cla$handangle_deg)
  
  attach(cla)
  
  if (factors == 'both')  {
    cat('\nLME with hand angle and movement type as fixed effects, and participant as random effect:\n\n')
    if (LMEmethod=='chi-squared') {
    print(Anova(lme(taperror_deg ~ passive_b * handangle_deg, random = ~1|participant, na.action=na.exclude), type=3))
    }
    if (LMEmethod=='Satterthwaite') {
      cla_model_lmer <- lmer(taperror_deg ~ passive_b * handangle_deg - (1|participant), na.action=na.exclude)
      print(anova(cla_model_lmer,ddf='Satterthwaite',type=3))
    }
  }
  if (factors == 'movementtype') {
    cat('\nLME with movement type as fixed effects - ignoring hand angle, and participant as random effect:\n\n')
    if (LMEmethod=='chi-squared') {
    print(Anova(lme(taperror_deg ~ passive_b, random = ~1|participant, na.action=na.exclude), type=3))
    }
    if (LMEmethod=='Satterthwaite') {
      cla_model_lmer <- lmer(taperror_deg ~ passive_b - (1|participant), na.action=na.exclude)
      print(anova(cla_model_lmer,ddf='Satterthwaite',type=3))
    }
  }
  if (factors == 'handangle') {
    cat('\nLME with movement type as fixed effects - ignoring hand angle, and participant as random effect:\n\n')
    if (LMEmethod=='chi-squared') {
    print(Anova(lme(taperror_deg ~ handangle_deg, random = ~1|participant, na.action=na.exclude), type=3))
    }
    if (LMEmethod=='Satterthwaite') {
      cla_model_lmer <- lmer(taperror_deg ~ handangle_deg - (1|participant), na.action=na.exclude)
      print(anova(cla_model_lmer,ddf='Satterthwaite',type=3))
    }
  }
  
  
  detach(cla)
  
  options('contrasts' <- default.contrasts)
  
}

# boxPlotLocalization <- function() {
#   
#   # par(mfrow=c(4,1),mar=c(2,3,1,1))
#   
#   for (group in c('exposure','classic','online')) {
#     
#     df <- getPointLocalization(group, difference=FALSE, verbose=FALSE)
#     
#     # print(str(df))
#     
#     for (rotated in c(0,1)) {
#       
#       for (passive in c(0,1)) {
#         
#         # subdf <- df[which(df$rotated_b == rotated & df$passive_b == passive),]
#         
#         # boxplot(taperror_deg ~ handangle_deg, data=df, axes=F, bty='n')
#         
#         for (target in c(15,25,35,45,55,65,75)) {
#           
#           subdf <- df[which(df$rotated_b == rotated & df$passive_b == passive & df$handangle_deg == target),]
#           
#           ppno <- subdf$participant
#           taperror <- subdf$taperror_deg
#           
#           idx <- which(abs(taperror - mean(taperror)) > (3 * sd(taperror)))
#           
#           if (length(idx) > 0) {
#             
#             cat(sprintf('%s %s %s %d deg\n', group, c('aligned','rotated')[rotated+1], c('active','passive')[passive+1], target))
#             print(ppno[idx])
#             
#           }
#           
#         }
#         
#       }
#       
#     }
#     
#   }
#   
# }

getPeakLocConfidenceInterval <- function(group, CIs=c(.95), movementtype='active', LRpart='all', selectPerformance=TRUE, iterations=100000, remove15=FALSE) {
  
  filename <- sprintf('LOC_peakCI_%s.csv', group)
  
  if (file.exists(filename)) {
    
    cat(sprintf('\nloading peak LOC generalization from file for: %s\n',toupper(group)))
    
    df <- read.csv(filename, stringsAsFactors=FALSE)
    
  } else {
    
    cat(sprintf('\nbootstrapping peak LOC generalization for: %s\n',toupper(group)))
    
    loc <- getPointLocalization(group, difference=TRUE, points=c(15,25,35,45,55,65,75), movementtype=movementtype, LRpart=LRpart, verbose=FALSE, selectPerformance=selectPerformance)
    
    if (remove15) {
      loc <- loc[which(loc$handangle_deg > 15),]
    }
    
    loc2 <- -1 * xtabs(taperror_deg ~ participant + handangle_deg,loc)
    
    data <- bootstrapGaussianPeak(data=loc2,bootstraps=iterations,mu=47.5,sigma=15,scale=10,offset=4,CIs=CIs)
    
    df <- data.frame('level'=names(data),'value'=data)
    
    write.csv(df,filename,row.names=FALSE,quote=FALSE)
    
  }
  
  return(df)
  
}

# DATA DESCRIPTIVES ------

countLocNAs <- function(group='exposure', output='count', selectPerformance=selectPerformance) {
  
  loc <- getPointLocalization(group, difference=FALSE, verbose=FALSE, selectPerformance=selectPerformance)
  loc <- loc[is.finite(loc$taperror_deg),]
  
  df <- expand.grid(unique(loc$participant), unique(loc$handangle_deg))
  names(df) <- c('participant', 'handangle_deg')
  df$count <- 0
  
  for (rown in c(1:nrow(df))) {
    
    pp    <- df[rown, 'participant']
    angle <- df[rown, 'handangle_deg']
    
    subexp <- loc[which(loc$participant == pp & loc$handangle_deg == angle),]
    df$count[rown] <- nrow(subexp)
    
  }
  
  if (output == 'count') {
    # this is a count of participants with estimates in ALL 4 tasks
    df$count[which(df$count < 4)] <- 0
    df$count[which(df$count > 0)] <- 1
    df <- aggregate(count ~ handangle_deg, data=df, FUN=sum)
  }
  
  if (output == 'percentage') {
    # this returns a percentage of existing estimates across the 4 tasks
    df$count <- df$count/4
    df <- aggregate(count ~ handangle_deg, data=df, FUN=mean)
  }
  
  return(df)
  
}

getLocCountTable <- function(output='count', selectPerformance=selectPerformance) {
  
  groups <- c('exposure','classic','online')
  
  for (group in groups) {
    
    counts <- countLocNAs(group=group, output=output, selectPerformance=selectPerformance)
    
    if (group == groups[1]) {
      
      df <- counts
      names(df)[2] <- group
      
    } else {
      
      df[,group] <- counts$count
      
    }
    
  }
  
  return(df)
  
}

countSelectedLocalizations <- function(group, ignoreRepetitions=FALSE, selectPerformance=TRUE) {
  
  df <- load.DownloadDataframe(url=localizationURLs[group],filename=sprintf('localization_%s.csv',group))
  
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
  
  mintrials <- 25
  
  participants <- unique(df$participant)
  
  for (ppid in participants) {
    
    ppdf <- df[which(df$participant == ppid),]
    
    for (session in c(0,1)) {
      
      for (movtype in c(0,1)) {
        
        subdf <- ppdf[which(ppdf$rotated_b == session & ppdf$passive_b == movtype),]
        
        iters <- unique(subdf$iteration)
        
        for (iterno in c(1:length(iters))) {
          
          iter <- iters[iterno]
          
          iterdf <- subdf[which(subdf$iteration == iter),]
          
          Ntrials <- dim(iterdf)[1]
          
          if (Ntrials < mintrials) {
            mintrials <- Ntrials
          }
          
          participant <- c(participant, ppid)
          rotated     <- c(rotated, session)
          passive     <- c(passive, movtype)
          repetition  <- c(repetition, iter)
          trials      <- c(trials, (Ntrials/.25))
          
        }
        
      }
      
    }
    
  }
  
  #cat(sprintf('\nminimum trials selected: %d\n\n',mintrials))
  
  return(data.frame(participant, rotated, passive, repetition, trials))
  
}