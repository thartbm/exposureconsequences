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

groupURLs <- c('exposure'='https://osf.io/9qfhp/download', 'classic'='https://osf.io/upw49/download', 'online'='https://osf.io/wjcgk/download')

# no-cursor reach data is not uploaded yet (not sure about the format)

localizationLMEsChiSq <- function() {
  
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
  

}

localizationLMEsSatterthwaite <- function() {
  
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
  
}

getANOVAlocalization <- function(group) {
  
  df <- read.csv(groupURLs[group],stringsAsFactors=FALSE)
  
  df <- aspligned(df)
  
  groupdf <- NA
  
  participants <- unique(df$participant)
  
  for (pp.no in c(1:length(participants))) {
    
    pp.id <- participants[pp.no]
    
    for (rotated in c(0,1)) {
      
      for (passive in c(0,1)) {
        
        subdf <- df[which(df$participant == pp.id & df$rotated_b == rotated & df$passive_b == passive),]
        
        locdf <- getLocalizationPoints(subdf, points=c(15,25,35,45,55,65,75), removeOutliers=TRUE)
        
        if (any(is.na(locdf$taperror_deg))) {
          cat(sprintf('WARNING: NAs in %s, pp:%d, %s/%s\n',group,pp.id,c('aligned','rotated')[rotated+1],c('active','passive')[passive+1]))
        }
        
        if (is.data.frame(groupdf)) {
          groupdf <- rbind(groupdf, locdf)
        } else {
          groupdf <- locdf
        }
        
      }
      
    }
    
  }
  
  return(groupdf)
  
}

aspligned <- function(df) {
  
  participants <- unique(df$participant)
  
  for (pp.no in c(1:length(participants))) {
    
    # fit aligned data with a smooth spline:
    al.idx <- which(df$participant == pp.no & df$rotated_b == 0)
    X <- df$handangle_deg[al.idx]
    Y <- df$taperror_deg[al.idx]
    AL.spl <- smooth.spline(x=X,y=Y,spar=0.99, keep.data=FALSE)
    
    # subtract predicted errors from both aligned and rotated data, using the fitted spline
    pp.idx <- which(df$participant == pp.no)
    predicted_errors <- predict(AL.spl, x=df$handangle_deg[pp.idx])
    df$taperror_deg[pp.idx] <- df$taperror_deg[pp.idx] - predicted_errors$y

  }
  
  return(df)
  
}

getLocalizationPoints <- function(subdf, points=c(15,25,35,45,55,65,75), removeOutliers=FALSE) {
  
  X <- subdf$handangle_deg
  Y <- subdf$taperror_deg
  
  # remove data points out of range:
  idx <- which(X > 0 & X < 90)
  X <- X[idx]
  Y <- Y[idx]
  
  # see if samples are reasonably close to what a smoothed spline on the rest of the data would predict:
  if (removeOutliers) {
    
    deviations <- c()
    
    for (idx.idx in c(1:length(X))) {
      
      sampleX <- X[idx.idx]
      sampleY <- Y[idx.idx]
      
      spl <- smooth.spline(X[-idx.idx],Y[-idx.idx], spar=.99, keep.data=FALSE) # spar=.95
      
      sampleP <- predict(spl, x=sampleX)
      deviations <- c(deviations, sampleP$y - sampleY)
      
    }
    
    idx <- which(abs(deviations - mean(deviations)) < (3 * sd(deviations)))
    
    X <- X[idx]
    Y <- Y[idx]
    
  }
  
  # cubic spline object, based on the data:
  spl <- smooth.spline(X, Y, spar=0.65, keep.data=FALSE) # spar=0.65
  
  # predict (interpolate) at given coordinates:
  spl.pr <- predict(spl,points)
  # spl.pr$y now has the estimated values of Y (the localization error) at the points of interest
  
  group <- rep(subdf$group[1],length(points))
  online_b <- rep(subdf$online_b[1],length(points))
  participant <- rep(subdf$participant[1],length(points))
  rotated_b <- rep(subdf$rotated_b[1],length(points))
  passive_b <- rep(subdf$passive_b[1],length(points))
  handangle_deg <- points
  taperror_deg <- round(spl.pr$y, digits=5)
  
  # since splines do well at interpolating, but not extrapolating
  # everything outside of the range of the data should be set to NAs:
  Xrange <- range(X)
  taperror_deg[which(points < Xrange[1])] <- NA
  taperror_deg[which(points > Xrange[2])] <- NA
  
  return(data.frame(group, online_b, participant, rotated_b, passive_b, handangle_deg, taperror_deg))
  
}
