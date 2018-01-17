
source('shared.R')



correlateNoCursorsLocalization <- function() {
  
  # first we do all responses, plan to split it by target later
  # get all data, for classic and exposure, we need localization and reach aftereffects
  
  # LOCALIZATION SHIFT for exposure
  local.exp <- getANOVAlocalization(group='exposure')
  loc.exp.ro <- aggregate(taperror_deg ~ participant + passive_b + handangle_deg, data=local.exp[local.exp$rotated_b == 1 & local.exp$passive_b == 0,], FUN=mean, drop=FALSE)
  loc.exp.al <- aggregate(taperror_deg ~ participant + passive_b + handangle_deg, data=local.exp[local.exp$rotated_b == 0 & local.exp$passive_b == 0,], FUN=mean, drop=FALSE)
  local.exp <- loc.exp.ro
  local.exp$taperror_deg <- local.exp$taperror_deg - loc.exp.al$taperror_deg
  
  # REACH AFTEREFFECTS for exposure
  nocur.exp <- getReachAftereffects(group='exposure', part='all', clean=TRUE) 
  names(nocur.exp)[which(names(nocur.exp) == 'target')] <- 'handangle_deg'
  
  # combine data frames for exposure
  exposure <- cbind(nocur.exp[which(nocur.exp$participant %in% local.exp$participant),], local.exp)
  
  # LOCALIZATION SHIFT for classic
  local.cla <- getANOVAlocalization(group='classic')
  loc.cla.ro <- aggregate(taperror_deg ~ participant + passive_b + handangle_deg, data=local.cla[local.cla$rotated_b == 1 & local.cla$passive_b == 0,], FUN=mean, drop=FALSE)
  loc.cla.al <- aggregate(taperror_deg ~ participant + passive_b + handangle_deg, data=local.cla[local.cla$rotated_b == 0 & local.cla$passive_b == 0,], FUN=mean, drop=FALSE)
  local.cla <- loc.cla.ro
  local.cla$taperror_deg <- local.cla$taperror_deg - loc.cla.al$taperror_deg
  
  # REACH AFTEREFFECTS for classic
  nocur.cla <- getReachAftereffects(group='classic', part='all', clean=TRUE) 
  names(nocur.cla)[which(names(nocur.cla) == 'target')] <- 'handangle_deg'
  
  # cbind the dataframes fro classic:
  classic <- cbind(local.cla, nocur.cla)
  
  # now do the tests:
  
  classic.cor.test <- cor.test(classic$endpoint_angle, classic$taperror_deg)
  print(classic.cor.test)
  
  exposure.cor.test <- cor.test(exposure$endpoint_angle, exposure$taperror_deg)
  print(exposure.cor.test)
  
  # make a figure:
  par(mfrow=c(1,2))
  
  plot(classic$endpoint_angle, classic$taperror_deg, main='classic', xlab='reach aftereffect', ylab='localization shift')
  # fit regression model
  LOC <- classic$taperror_deg
  RAE <- classic$endpoint_angle
  classic.lm <- lm(LOC ~ RAE)
  # and use that to show a line:
  abline(classic.lm$coefficients, col='red')
  pointlocs <- seq(min(RAE),max(RAE),.1)
  
  lines(
    x   = pointlocs,
    y   = predict( classic.lm, 
                   newdata=data.frame(RAE=pointlocs), 
                   interval = "confidence" )[ , "upr" ],
    col = "red",
    lty = 2)
  
  lines(
    x   = pointlocs,
    y   = predict( classic.lm, 
                   newdata=data.frame(RAE=pointlocs), 
                   interval = "confidence" )[ , "lwr" ],
    col = "red",
    lty = 2)
  
  
  plot(exposure$endpoint_angle, exposure$taperror_deg, main='exposure', xlab='reach aftereffect', ylab='localization shift')
  # fit regression model
  LOC <- exposure$taperror_deg
  RAE <- exposure$endpoint_angle
  exposure.lm <- lm(LOC ~ RAE)
  # and use that to show a line:
  abline(exposure.lm$coefficients, col='blue')
  pointlocs <- seq(min(RAE),max(RAE),.1)
  
  lines(
    x   = pointlocs,
    y   = predict( exposure.lm, 
                   newdata=data.frame(RAE=pointlocs), 
                   interval = "confidence" )[ , "upr" ],
    col = "blue",
    lty = 2)
  
  lines(
    x   = pointlocs,
    y   = predict( exposure.lm, 
                   newdata=data.frame(RAE=pointlocs), 
                   interval = "confidence" )[ , "lwr" ],
    col = "blue",
    lty = 2)
  
  
}