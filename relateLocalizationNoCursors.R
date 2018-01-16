
source('shared.R')



correlateNoCursorsLocalization <- function() {
  
  # first we do all responses, then split by target
  
  # get all data, for classic and exposure, we need localization and reach aftereffects
  
  local.exp <- getANOVAlocalization(group='exposure')
  nocur.exp <- getReachAftereffects(group='exposure', part='all', clean=TRUE) 
  
  local.cla <- getANOVAlocalization(group='classic')
  nocur.cla <- getReachAftereffects(group='classic', part='all', clean=TRUE) 
  
  
  
}