
source('shared.R')

# first make sure we have all necessary packages installed and loaded:
# for now we'll use Chi-square (not Satterthwaite estimates of F/p-values):
installRequire.Packages(c('nlme', 'car'))

# , 'lme4', 'lmerTest'

# localization data can be downloaded from OSF:
# groupURLs <- c('exposure'='https://osf.io/9qfhp/download', 'classic'='https://osf.io/upw49/download', 'online'='https://osf.io/wjcgk/download')

# no-cursor reach data is not uploaded yet (not sure about the format)

localizationLMEsChiSq <- function() {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))

  exposure <- getANOVAlocalization('exposure')
  classic <- getANOVAlocalization('classic')
  
  exposure$participant <- exposure$participant + 100

  exposure <- exposure[-which(exposure$handangle_deg == 15),]
  classic <- classic[-which(classic$handangle_deg == 15),]
  
  # they all have this header:
  # group,online_b,participant,rotated_b,passive_b,handangle_deg,tapangle_deg
  # (where group is uninformative)

  localization <- rbind(exposure,classic)
  
  localization$group         <- factor(localization$group)
  localization$participant   <- factor(localization$participant)
  localization$online_b      <- factor(localization$online_b)
  localization$rotated_b     <- factor(localization$rotated_b)
  localization$passive_b     <- factor(localization$passive_b)
  localization$handangle_deg <- factor(localization$handangle_deg)
  
  # str(localization)
  
  # prepare a data frame with the training-induced changes in hand localization:
  locshift <- localization[which(localization$rotated_b == 1),]
  locshift$taperror_deg <- locshift$taperror_deg - localization$taperror_deg[which(localization$rotated_b == 0)]
  
  # first a model including the aligned and rotated data
  # here we look for any effect of session (rotated_b)
  # and if it is there, we continue with other analyses that are easier to grasp:
  exp_full_model <- lme(taperror_deg ~ rotated_b * passive_b * group, data=localization, random = ~1|participant/handangle_deg, na.action=na.exclude)
  print(Anova(exp_full_model, type=3))

  # the main question is whether or not the training type ('group') interacts with:
  # 1) the effect of the type of movement, active or passive ('passive_b') before localization
  # and perhaps on:
  # 2) the location of the hand in the workspace ('handangle_deg')
  # so we make the simplest model that allows testing those two interactions:
  simple_exp_model <- lme(taperror_deg ~ group + handangle_deg:group + passive_b:group, data=locshift, random = ~1|participant, na.action=na.exclude)
  print(Anova(simple_exp_model, type=3))
  
  # given that one or two of the interactions in the simple model were significant,
  # we now want to test the main effects of the factors of interest in each group
  # (as a follow-up / post-hoc if you will)
  exposure_hand_model <- lme(taperror_deg ~ handangle_deg + passive_b, data=exposure, random = ~1|participant, na.action=na.exclude)
  print(Anova(exposure_hand_model, type=3))
  classic_hand_model <- lme(taperror_deg ~ handangle_deg + passive_b, data=classic, random = ~1|participant, na.action=na.exclude)
  print(Anova(classic_hand_model, type=3))
  
  options('contrasts' <- default.contrasts)
  
}

localizationLMEsSatterthwaite <- function() {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  exposure <- getANOVAlocalization('exposure')
  classic <- getANOVAlocalization('classic')
  
  exposure$participant <- exposure$participant + 100
  
  exposure <- exposure[-which(exposure$handangle_deg == 15),]
  classic <- classic[-which(classic$handangle_deg == 15),]
  
  # they all have this header:
  # group,online_b,participant,rotated_b,passive_b,handangle_deg,tapangle_deg
  # (where group is uninformative)
  
  localization <- rbind(exposure,classic)
  
  localization$group         <- factor(localization$group)
  localization$participant   <- factor(localization$participant)
  localization$online_b      <- factor(localization$online_b)
  localization$rotated_b     <- factor(localization$rotated_b)
  localization$passive_b     <- factor(localization$passive_b)
  localization$handangle_deg <- factor(localization$handangle_deg)
  
  # str(localization)
  
  locshift <- localization[which(localization$rotated_b == 1),]
  locshift$taperror_deg <- locshift$taperror_deg - localization$taperror_deg[which(localization$rotated_b == 0)]
  
  # "omnibus" test
  exp_full_model_lmer <- lmer(taperror_deg ~ rotated_b * passive_b * group - (1|participant/handangle_deg), data=localization, na.action=na.exclude)
  print(anova(exp_full_model_lmer,ddf='Satterthwaite',type=3))
  
  # with restricted model:
  simple_exp_model_lmer <- lmer(taperror_deg ~ group+ handangle_deg:group + passive_b:group - (1|participant), data=locshift, na.action=na.exclude)
  print(anova(simple_exp_model_lmer,ddf='Satterthwaite',type=3))
  
  # split by group:
  exposure_hand_model_lmer <- lmer(taperror_deg ~ handangle_deg + passive_b - (1|participant), data=exposure, na.action=na.exclude)
  print(anova(exposure_hand_model_lmer,ddf='Satterthwaite',type=3))
  classic_hand_model_lmer <- lmer(taperror_deg ~ handangle_deg + passive_b - (1|participant), data=classic, na.action=na.exclude)
  print(anova(classic_hand_model_lmer,ddf='Satterthwaite',type=3))
  
  options('contrasts' <- default.contrasts)
  
}




# PLOTS / FIGURES ------

plotLocalization <- function() {
  
  # get the data to plot:
  exp <- getPointLocalization('exposure', difference=TRUE, verbose=FALSE)
  cla <- getPointLocalization('classic', difference=TRUE, verbose=FALSE)
  
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

  
  par(mfrow=c(1,2))
  points <- c(15,25,35,45,55,65,75)
  
  # panel A: exposure localization (active vs. passive)
  
  plot(-1000,-1000, main='exposure', xlab='hand angle [deg]', ylab='localization shift [deg]', xlim=c(10,80), ylim=c(0,-15), axes=F)
  
  X <- c(points, rev(points))
  exp.act.Y <- c(exp.CI.act[1,],rev(exp.CI.act[2,]))
  exp.pas.Y <- c(exp.CI.pas[1,],rev(exp.CI.pas[2,]))
  
  polygon(X,exp.act.Y,border=NA,col=colorset[['expActT']])
  polygon(X,exp.pas.Y,border=NA,col=colorset[['expPasT']])
  
  lines(points[1:2],exp.avg.act$taperror_deg[1:2],col=colorset[['expActS']],lty=2,lwd=2)
  lines(points[2:7],exp.avg.act$taperror_deg[2:7],col=colorset[['expActS']],lty=1,lwd=2)
  lines(points[1:2],exp.avg.pas$taperror_deg[1:2],col=colorset[['expPasS']],lty=2,lwd=2)
  lines(points[2:7],exp.avg.pas$taperror_deg[2:7],col=colorset[['expPasS']],lty=1,lwd=2)
  
  axis(1,at=points)
  axis(2,at=c(0,-5,-10,-15))
  
  legend(10,-15,c('passive','active'),col=c(colorset[['expPasS']],colorset[['expActS']]),lty=c(1,1),lwd=c(2,2),bty='n')
  
  # panel B: classic localization (active vs. passive)
  
  plot(-1000,-1000, main='classic', xlab='hand angle [deg]', ylab='localization shift [deg]', xlim=c(10,80), ylim=c(0,-15), axes=F)
  
  X <- c(points, rev(points))
  cla.act.Y <- c(cla.CI.act[1,],rev(cla.CI.act[2,]))
  cla.pas.Y <- c(cla.CI.pas[1,],rev(cla.CI.pas[2,]))
  
  polygon(X,cla.act.Y,border=NA,col=colorset[['claActT']])
  polygon(X,cla.pas.Y,border=NA,col=colorset[['claPasT']])
  
  lines(points[1:2],cla.avg.act$taperror_deg[1:2],col=colorset[['claActS']],lty=2,lwd=2)
  lines(points[2:7],cla.avg.act$taperror_deg[2:7],col=colorset[['claActS']],lty=1,lwd=2)
  lines(points[1:2],cla.avg.pas$taperror_deg[1:2],col=colorset[['claPasS']],lty=2,lwd=2)
  lines(points[2:7],cla.avg.pas$taperror_deg[2:7],col=colorset[['claPasS']],lty=1,lwd=2)
  
  axis(1,at=points)
  axis(2,at=c(0,-5,-10,-15))
  
  legend(10,-15,c('passive','active'),col=c(colorset[['claPasS']],colorset[['claActS']]),lty=c(1,1),lwd=c(2,2),bty='n')
  
  
}

exposureLocalization <- function(remove15=TRUE) {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  exp <- getPointLocalization('exposure', difference=FALSE, verbose=FALSE)
  
  if (remove15) {
    exp <- exp[-which(exp$handangle_deg == 15),]
  }
  
  exp$participant   <- factor(exp$participant)
  exp$rotated_b    <- factor(exp$rotated_b)
  exp$passive_b     <- factor(exp$passive_b)
  exp$handangle_deg <- factor(exp$handangle_deg)
  
  attach(exp)
  
  cat('\nLME with session, target and movement type as fixed effects, and participant as random effect:\n\n')
  print(Anova(lme(taperror_deg ~ rotated_b * passive_b * handangle_deg, random = ~1|participant, na.action=na.exclude), type=3))
  
  detach(exp)
  
  options('contrasts' <- default.contrasts)
  
}

exposureLocalizationShift <- function(noHandAngle=FALSE, remove15=TRUE) {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  exp <- getPointLocalization('exposure', difference=TRUE, verbose=FALSE)
  
  if (remove15) {
    exp <- exp[-which(exp$handangle_deg == 15),]
  }
  
  exp$participant   <- factor(exp$participant)
  exp$passive_b     <- factor(exp$passive_b)
  exp$handangle_deg <- factor(exp$handangle_deg)
  
  attach(exp)
  
  if (noHandAngle) {
    cat('\nLME with movement type as fixed effects - ignoring hand angle, and participant as random effect:\n\n')
    print(Anova(lme(taperror_deg ~ passive_b, random = ~1|participant, na.action=na.exclude), type=3))
    
  } else {
    cat('\nLME with hand angle and movement type as fixed effects, and participant as random effect:\n\n')
    print(Anova(lme(taperror_deg ~ passive_b * handangle_deg, random = ~1|participant, na.action=na.exclude), type=3))
  }
  
  detach(exp)
  
  options('contrasts' <- default.contrasts)
  
}

groupLocalization <- function(model='full',remove15=TRUE) {
  
  default.contrasts <- options('contrasts')
  options(contrasts=c('contr.sum','contr.poly'))
  
  exp <- getPointLocalization('exposure', difference=TRUE, verbose=FALSE)
  cla <- getPointLocalization('classic', difference=TRUE, verbose=FALSE)
  
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
    print(Anova(lme(taperror_deg ~ group * passive_b * handangle_deg, random = ~1|participant, na.action=na.exclude), type=3))
  }
  if (model == 'restricted') {
    cat('\nLME with three terms only, removing some main effects and interactions:\n\n')
    print(Anova(lme(taperror_deg ~ group + group:passive_b + group:handangle_deg, random = c(~1|passive_b, ~1|handangle_deg, ~1|participant), na.action=na.exclude), type=3))
  }
  if (model == 'handangle') {
    cat('\nLME with group and hand angle as fixed effects, and participant and movement type as random effects:\n\n')
    print(Anova(lme(taperror_deg ~ group * handangle_deg, random = ~1|participant/taperror_deg, na.action=na.exclude), type=3))
  }
  if (model == 'movementtype') {
    cat('\nLME with group and movement type as fixed effects and participant and hand angle as random effects:\n\n')
    print(Anova(lme(taperror_deg ~ group * passive_b, random = ~1|participant/handangle_deg, na.action=na.exclude), type=3))
  }
  
  detach(loc)
  
  options('contrasts' <- default.contrasts)
  
}

classicLocalizationShift <- function(factors='both',remove15=TRUE) {
  
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
    print(Anova(lme(taperror_deg ~ passive_b * handangle_deg, random = ~1|participant, na.action=na.exclude), type=3))
  }
  if (factors == 'movementtype') {
    cat('\nLME with movement type as fixed effects - ignoring hand angle, and participant as random effect:\n\n')
    print(Anova(lme(taperror_deg ~ passive_b, random = ~1|participant, na.action=na.exclude), type=3))
  }
  if (factors == 'handangle') {
    cat('\nLME with movement type as fixed effects - ignoring hand angle, and participant as random effect:\n\n')
    print(Anova(lme(taperror_deg ~ handangle_deg, random = ~1|participant, na.action=na.exclude), type=3))
  }
  
  
  detach(cla)
  
  options('contrasts' <- default.contrasts)
  
}

boxPlotLocalization <- function() {
  
  # par(mfrow=c(4,1),mar=c(2,3,1,1))
  
  for (group in c('exposure','classic','online')) {
    
    df <- getPointLocalization(group, difference=FALSE, verbose=FALSE)
    
    # print(str(df))
    
    for (rotated in c(0,1)) {
      
      for (passive in c(0,1)) {
        
        # subdf <- df[which(df$rotated_b == rotated & df$passive_b == passive),]
        
        # boxplot(taperror_deg ~ handangle_deg, data=df, axes=F, bty='n')
        
        for (target in c(15,25,35,45,55,65,75)) {
          
          subdf <- df[which(df$rotated_b == rotated & df$passive_b == passive & df$handangle_deg == target),]
          
          ppno <- subdf$participant
          taperror <- subdf$taperror_deg
          
          idx <- which(abs(taperror - mean(taperror)) > (3 * sd(taperror)))
          
          if (length(idx) > 0) {
            
            cat(sprintf('%s %s %s %d deg\n', group, c('aligned','rotated')[rotated+1], c('active','passive')[passive+1], target))
            print(ppno[idx])
            
          }
          
        }
        
      }
      
    }
    
  }
  
}