######################################################################
##### Code for Day 2 of SEM Workshop
##### covering assessing SEMs of observed variables
#####
##### Jarrett E.K. Byrnes
#####
##### Last Modified 10/13/17
######################################################################


#Load your Data File
keeley<-read.csv("../data/Keeley_rawdata_select4.csv")

#load lavaan
library(lavaan)

#load library for multivarite normality test of residuals
#library(mvnormtest)

#load some helper functions for along the way
#source("./fitted_lavaan.R")

####################################
# Evaluating Model Fit          ####
####################################


fullMedModel<-' firesev ~ age
cover ~ firesev'

fullMedSEM<-sem(fullMedModel, 
                data=keeley, 
                meanstructure=TRUE)

#See the model chi-square
fullMedSEM

#get fit statistics.  All of 'em
summary(fullMedSEM, fit.measures=T)
fitMeasures(fullMedSEM)

#######################################
# What did you miss in your model? ####
#######################################

#Sample versus model covariance
inspect(fullMedSEM, "sample")
inspect(fullMedSEM, "fitted")


#Residuals
residuals(fullMedSEM)

#correlation matrix
residuals(fullMedSEM, type="cor")


#Modification Indices
modificationIndices(fullMedSEM, standardized=F)

##########################################
# Exercise: Diagnose misspecification ####
##########################################


#The Full Mediation Model
distMedModel <- '
  rich ~ abiotic + hetero
  hetero ~ distance
  abiotic ~ distance'

distMedFit <- sem(distMedModel, data=keeley)

distMedFit

#Look at residuals
residuals(distMedFit, type="cor")

#Look at Lagrange multipliers
modI<-modificationIndices(distMedFit, standardized=F)
modI[modI$mi>3,]



####################################
# Evaluating Residuals
####################################
source("./fitted_lavaan.R")

#The Richness Partial Mediation Model
distModel <- 'rich ~ distance + abiotic + hetero
hetero ~ distance
abiotic ~ distance'

distFit <- sem(distModel, data=keeley)

res <- residuals(distFit, "casewise")

distFitInt <- sem(distModel, data=keeley, meanstructure=T)

res <- residuals_lavaan(distFitInt)

head(res)

par(mfrow=c(1,3))

apply(res[,1:3], 2, function(x){
  qqnorm(x, cex=1.5, cex.lab=1.5, cex.axis=1.3)
  qqline(x)
})

par(mfrow=c(1,1))


#####################################
# Testing Normality of residuals ####
#####################################
#Some code to help you out
#will be integrated into semTools forthwith!
source("./fitted_lavaan.R")

#Fit model
partialMedModel<-' firesev ~ age
                cover ~ firesev + age'

partialMedSEM<-sem(partialMedModel, 
                   data=keeley,
                   meanstructure=T)

#get residuals
partialResid <- residuals_lavaan(partialMedSEM)

#see those residuals
head(partialResid)


#Plot all the residuals
qqnorm_plot <- function(x){qqnorm(x); qqline(x)}

#Two panels
par(mfrow=c(1,2))
apply(partialResid, 2, qqnorm_plot)
par(mfrow=c(1,1))

#Test with Shapiro-Wilks
library(mvnormtest)

mshapiro.test(t(partialResid))

#Or, Mardia's test of MV Normality
library(MVN)


mvn(partialResid, mvnTest="mardia")

mvn(partialResid, mvnTest="mardia",  univariatePlot = "qqplot")

mvn(partialResid, mvnTest="mardia",  multivariatePlot = "persp")


#########################
# Residuals Exercise ####
#########################

distMedFit <- sem(distMedModel,
                  data=keeley,
                  meanstructure=TRUE)

dist_resid <- residuals_lavaan(distMedFit)

#plot it and evaluate
mvn(dist_resid, mvnTest="mardia",  univariatePlot = "qqplot")


####################################
# Non-normality of data         ####
####################################

library(scatterplot3d)
fitdata <- inspect(distMedFit, "data")

scatterplot3d(fitdata)

#test for multivariate normality
mvn(fitdata, mvnTest="mardia",  univariatePlot = "qqplot")


#################################################
# Adjusting to Non-normality of data         ####
#################################################

#SB Chisq
distFitSB<-sem(distModel, data=keeley, estimator="MLM")
distFitSB

#Bootstrap
distFitBoot<-sem(distModel, data=keeley, 
                  test="bollen.stine", se="boot", bootstrap=100)
distFitBoot


####################################
# Testing Mediation             ####
####################################

#The models
fullMedModel<-' firesev ~ age
                cover ~ firesev'

fullMedSEM<-sem(fullMedModel,
                data=keeley)

partialMedModel<-' firesev ~ age
                cover ~ firesev + age'

partialMedSEM<-sem(partialMedModel, 
                   data=keeley)


#with a LRT
anova(partialMedSEM, fullMedSEM)

#with aicc
library(AICcmodavg)


aictab(cand.set = list(fullMedSEM, partialMedSEM), 
       modnames = c("Full", "Partial"))

###
#
###


#The  Partial Mediation Model
distModel <- 'rich ~ distance + abiotic + hetero
hetero ~ distance
abiotic ~ distance'

distFit <- sem(distModel, data=keeley)


#The Full Mediation Model
distMedModel <- '
rich ~ abiotic + hetero
hetero ~ distance
abiotic ~ distance'

distMedFit <- sem(distMedModel, data=keeley)

#LRT
anova(distFit, distMedFit)

#AIC
aictab(cand.set = list(distMedFit, distFit), 
       modnames = c("Full", "Partial"))


#############
#POWER
#############
library(simsem)

buildPopMod <- function(fit){
  ptab <- parTable(fit) 
  if(length(which(ptab$exo==1))>0)
    ptab <- ptab[-which(ptab$exo==1),]
  mod <- with(ptab, paste(lhs, op, est, "*", rhs, collapse="\n\n"))
  mod
}

n <- round(seq(5,1000, length.out=100))
p <- rep(NA, length(n))
for(i in 1:length(n))
  p[i] <- SSpower(buildPopMod(distMedFit), n[i], distModel, fun="sem")

plot(n, p)



