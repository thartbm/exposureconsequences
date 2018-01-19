
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

