##################
# Code for analysis of categorical data
# in SEMs
#
# Jarrett Byrnes
##################

library(piecewiseSEM)
library(nlme)

#data from https://github.com/jebyrnes/phrag_common_garden_sem
#Bowen, J.L., Kearns, P.J., Byrnes, J.E.K., Wigginton, S., Allen, W.J., 
#Greenwood, M., Tran, K., Yu, J., Cronin, J.T., Meyerson, L.A., 2017. Lineage 
#overwhelms environmental conditions in determining rhizosphere bacterial 
#community structure in a cosmopolitan invasive plant. Nature Communications 
#8, 501. 

bowen <- read.csv("../data/bowen.csv")

###
# A simple categorical model
###
div_mod <- lme(observed_otus ~  status, 
                      random =~ 1|Genotype,
                      data = bowen, method = "ML")

activity_mod <- lme(RNA.DNA ~ status + observed_otus, 
                           random =~ 1|Genotype, 
                           data=bowen, method="ML")

c_mod <- lme(below.C ~ observed_otus +  status, 
                    random =~ 1|Genotype, data=bowen, method="ML")

biomass_mod <- lme(abovebiomass_g ~ RNA.DNA + observed_otus + below.C + status, 
                           random =~ 1|Genotype, data = bowen, method="ML")

bowen_mod <- psem(
  div_mod,
  activity_mod,
  c_mod,
  biomass_mod,
  data = bowen
)

basisSet(bowen_mod)
dSep(bowen_mod)
coefs(bowen_mod)

#what do the ANOVAs say?
anova(bowen_mod)

#The coefficients
library(emmeans)

lapply(bowen_mod[-length(bowen_mod)], emmeans, specs = ~status )


###
# Test of Differences of means
###

#Let's look at, posthoc tests
generic_tukey <- function(x)
  emmeans(x, list(pairwise ~ status))


lapply(bowen_mod[-length(bowen_mod)], generic_tukey)

####
## Multigroup ####
####

meadows<-read.csv("../data/FinnishMeadows_Multigroup.csv")


####################################################
# Multigroup Analysis with One Coefficient      ####
####################################################



meadowModel<-'rich ~ elev + mass
mass ~ me*elev'


#different fits - one whole shebang, once grouped with everything varying, one grouped, but with everything constrained
meadowFit<-sem(meadowModel, data=meadows)

coef(meadowFit)
summary(meadowFit)

#OK, now introduce groups
meadowFitFree<-sem(meadowModel, data=meadows, group="grazed")
coef(meadowFitFree)

meadowFitFree
summary(meadowFitFree)

#introducing constraints
meadowFitEqual<-sem(meadowModel, data=meadows, 
                    group="grazed", 
                    group.equal=c( "intercepts","regressions"))

coef(meadowFitEqual)

####
# Comparing constrained and unconstrained models ####
####


#is the constraint OK?
anova(meadowFitFree, meadowFitEqual)

#constrain just the elev-mass relationship
meadowModel2<-'rich ~ c("me", "me")*elev + mass
mass ~ elev'

meadowFitFree2<-sem(meadowModel2, data=meadows, group="grazed")

coef(meadowFitFree2)
summary(meadowFitFree2)

#test to see if this is one of the culprits!
anova(meadowFitFree, meadowFitFree2)

###############
#Exercise- split by par2 ####
###############

meadowFitFreePar2<-sem(meadowModel, data=meadows, group="par2")
meadowFitEqualPar2<-sem(meadowModel, data=meadows, group="par2", 
                        group.equal=c( "intercepts", "regressions"))

summary(meadowFitFreePar2)

# Test Constraints
anova(meadowFitEqualPar2, meadowFitFreePar2)


######
# Exercise: Try other constraints & group - here's one
######

meadowModel3<-'rich ~ c("re", "re")*elev + mass
mass ~ elev'

meadowFitFreePar2_3<-sem(meadowModel3, data=meadows, group="par2")

#does it matter?  Yes!
anova(meadowFitFreePar2, meadowFitFreePar2_3)

#see others that may vary
coef(meadowFitFreePar2_3)

######
# Exercise: Try releasing constraints
######


meadowModelNoLabel<-'rich ~ elev + mass
mass ~ elev'


meadowFitOneFree<-sem(meadowModelNoLabel, data=meadows, group="grazed", 
                      group.equal = c("regressions", "intercepts"),
                      group.partial = c("mass ~ elev", "mass ~ 1"))


summary(meadowFitOneFree)

meadowFitTwoFree<-sem(meadowModelNoLabel, data=meadows, group="grazed", 
                      group.equal = c("regressions", "intercepts"),
                      group.partial=c("mass ~ elev", "mass ~ 1",
                                      "rich ~ mass", "rich ~ 1"))

summary(meadowFitTwoFree)


#compare to constrained model
anova(meadowFitTwoFree, meadowFitFree)


###
# Multigroup with piecewiseSEM ###
###

meadows$grazed <- factor(meadows$grazed)

#Fully unconstrained model
rich_unconstrained <- lm(rich ~ elev*grazed + mass * grazed, data=meadows)
mass_unconstrained <- lm(mass ~ elev * grazed, data=meadows)

unconstrained_int_mod <- psem(
  rich_unconstrained,
  mass_unconstrained,
  data=meadows
)

anova(unconstrained_int_mod)


# One Constraint
rich_constrained <- lm(rich ~ elev + mass * grazed, data=meadows)

constrained_int_mod <- psem(
  rich_constrained,
  mass_unconstrained,
  data=meadows
)

anova(unconstrained_int_mod, constrained_int_mod)

source("./multigroup_margins.R")
getMargins.psem(unconstrained_int_mod, 
                at=list(grazed=c(0,1)))

getStdByGroup(unconstrained_int_mod, by = "grazed")
#####


