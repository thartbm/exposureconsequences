
source('shared.R')



correlateNoCursorsLocalization <- function() {
  
  # first we do all responses, plan to split it by target later
  # get all data, for classic and exposure, we need localization and reach aftereffects
  
  # LOCALIZATION SHIFT for exposure
  local.exp <- getANOVAlocalization(group='exposure')
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
  nocur.exp <- getReachAftereffects(group='exposure', part='all', clean=TRUE) 
  names(nocur.exp)[names(nocur.exp) == 'endpoint_angle'] <- 'RAE'
  # combine data frames for exposure
  exposure <- nocur.exp[which(nocur.exp$participant %in% local.exp.act$participant),]
  exposure$active_localization_shift <- local.exp.act$taperror_deg
  exposure$passive_localization_shift <- local.exp.pas$taperror_deg
  
  # LOCALIZATION SHIFT for classic
  local.cla <- getANOVAlocalization(group='classic')
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
  nocur.cla <- getReachAftereffects(group='classic', part='all', clean=TRUE) 
  names(nocur.cla)[names(nocur.cla) == 'endpoint_angle'] <- 'RAE'
  
  # cbind the dataframes fro classic:
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
  xlimits <- c(-5.5,25.5)
  ylimits <- c(-20.5,10.5)
  xticks <- seq(from=-5,to=25,by=5)
  yticks <- seq(from=-20,to=10,by=5)
  
  par(mfrow=c(2,2),mar=c(5.1,4.1,2.1,1.1))
  
  dataframes = list()
  dataframes[['exposure']] <- exposure
  dataframes[['classic']] <- classic
  
  for (group in c('exposure','classic')) {
    for (movement in c('active','passive')) {
      df <- dataframes[[group]]
      colname <- sprintf('%s_localization_shift',movement)
      X <- df$RAE
      Y <- df[,colname]
      
      this.cor.test <- cor.test(X, Y)
      
      cat(sprintf('%s %s\n',group,movement))
      print(this.cor.test)
      
      if (group == 'exposure') {
        colors = c(rgb(.5,.5,.5,.5),rgb(0,0,1))
      }
      if (group == 'classic') {
        colors = c(rgb(.5,.5,.5,.5),rgb(1,0,0))
      }
      main <- sprintf('%s %s',group,movement)
      xlab <- 'reach aftereffect'
      ylab <- sprintf('%s localization shift',movement)
      
      plotCorrelationWithCI(X=X,Y=Y,colors=colors,main=main,xlab=xlab,ylab=ylab,axes=FALSE,xticks=xticks,yticks=yticks,xlim=xlimits,ylim=ylimits)
      
    }
  }
  
}

plotCorrelationWithCI <- function(X,Y,colors=c('black','red'),main='main title',xlab='x',ylab='y',axes=FALSE,xticks=c(),yticks=c(),xlim=NULL,ylim=NULL) {
  
  plot(X, Y, col=colors[1], main=main, xlab=xlab, ylab=ylab, asp=1, axes=axes, xlim=xlim, ylim=ylim)
  # fit regression model
  
  lines(c(0,0),ylim,col=colors[1])
  lines(xlim,c(0,0),col=colors[1])
  
  
  this.lm <- lm(Y ~ X)
  # and use that to show a line:
  abline(this.lm$coefficients, col=colors[2])
  
  # now show the confidence interval
  pointlocs <- seq(min(X),max(X),.1)
  
  lines(
    x   = pointlocs,
    y   = predict( this.lm, 
                   newdata=data.frame(X=pointlocs), 
                   interval = "confidence" )[ , "upr" ],
    col = colors[2],
    lty = 2)
  
  lines(
    x   = pointlocs,
    y   = predict( this.lm, 
                   newdata=data.frame(X=pointlocs), 
                   interval = "confidence" )[ , "lwr" ],
    col = colors[2],
    lty = 2)
  
  if (length(xticks)) {
    axis(1,at=xticks)
  }
  if (length(yticks)) {
    axis(2,at=yticks)
  }
}