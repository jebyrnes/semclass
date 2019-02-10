######################################################################
##### Code for Day 3 of SEM Workshop
##### covering  Latent Variables
#####
##### Jarrett E.K. Byrnes
#####
##### Last Modified 2019-02-10
######################################################################

#Set your working Directory
library(lavaan)



####
#Multiple Indicators - Santos Example ####
####

santosCov <-read.table("../data/santosCov.txt", na.strings=".")
santosCov <-as.matrix(santosCov)

santosCFA1<-'
  Aposematism =~Alkaloid.quantity+Alkaloid.diversity+Conspicuous.coloration
'

santosFit1<-sem(santosCFA1, sample.cov=santosCov, sample.nobs=21)
summary(santosFit1)
standardizedSolution(santosFit1)	

####
# Exercise - body size ####
####

santosSize<-'
  Size =~ Log.Mass + Log.RMR + Log.Scope'

santosSizeFit<-sem(santosSize, sample.cov=santosCov, sample.nobs=21)

summary(santosSizeFit, standardized=T)

####
#Two variable CFA ####
####

santosCFA2<-paste(santosCFA1, '
		Aposematism =~ Ant.Mite.Specialization+log.Prey
                  Scale =~Log.Mass+Log.RMR+Log.Scope+Conspicuous.coloration
                  ', sep="\n")
	

santosFit<-sem(santosCFA2, sample.cov=santosCov, sample.nobs=21)
summary(santosFit)
standardizedSolution(santosFit)


####
#Spartina SEM ####
####
spartina <- read.csv("../data/TravisDataForLVExample.csv")


#what if an exogenous variable affects indicators of a latent variable
spartinaModel<-'performance =~  clonediam + numbstems + numbinfl + meanht + meanwidth
               meanht ~~ meanwidth

               performance ~ geneticdist
'

spartinaFit<-sem(spartinaModel, data=spartina)

summary(spartinaFit, standardized=T, rsquare=T)

############################
## Spartina Exercise ####
############################

spartinaModel2<-paste(spartinaModel, '
                      meanht ~ latitude
                      meanwidth ~ latitude', sep="\n")

spartinaFit2<-sem(spartinaModel2, data=spartina)

summary(spartinaFit2, standardized=T, rsquare=T)

############################
# Kelp Measurement Error ####
############################
library(ggplot2)

lter<-read.csv("../data/lter_kelp.csv")

#1) Calculate fitted values for spring biomass
#landsat observations to biomass
lter$landsat_spring_biomass<-154.89*lter$spring_canopy+68.62

#2) Calculate fitted values for summer biomass
#summer kelp counts to biomass y=0.08x+0.01 r^2=0.79
lter$summer_kelp_biomass<-0.08*lter$kelp+0.01


#3) Transform fitted values for easier fitting
#transformation for easier fitting
lter$summer_kelp_biomass<-log(lter$summer_kelp_biomass +1)
lter$landsat_spring_biomass <-log(lter$landsat_spring_biomass +1)

#No measurement error
noerror<-'summer_kelp_biomass ~ landsat_spring_biomass'
noerrorFit<-sem(noerror, data=lter)
summary(noerrorFit, standardized=T)

qplot(landsat_spring_biomass, summer_kelp_biomass, data=lter) +
  theme_bw(base_size=17)


#measurement error in predictor
var(lter$landsat_spring_biomass, na.rm=T)*(0.38)

#[1] 3.762301

errorCanopy<-'
  true_spring_biomass =~ 1*landsat_spring_biomass

  summer_kelp_biomass ~ true_spring_biomass

#error
landsat_spring_biomass ~~ 3.762301*landsat_spring_biomass
'

errorCanopyFit<-sem(errorCanopy, data=lter)
summary(errorCanopyFit, standardized=T)

error_pred <- as.data.frame(cbind(lavPredict(errorCanopyFit), inspect(errorCanopyFit, "data")))
names(error_pred) <- c("true_spring_biomass", "landsat_spring_biomass", "summer_kelp_biomass")

qplot(true_spring_biomass, summer_kelp_biomass, data=error_pred) +
  theme_bw(base_size=17)


##error in both
var(lter$summer_kelp_biomass, na.rm=T)*(0.21)
#[1] 0.5495345

errorBoth<-'
  true_spring_biomass =~ 1*landsat_spring_biomass
	true_summer_biomass =~ 1*summer_kelp_biomass
	
	true_summer_biomass ~ true_spring_biomass
	
	#error
	landsat_spring_biomass ~~ 3.762301*landsat_spring_biomass
	summer_kelp_biomass ~~ 0.5495345* summer_kelp_biomass
'

errorBothFit<-sem(errorBoth, data=lter)
summary(errorBothFit, standardized=T)


errorCanopyFit<-sem(errorCanopy, data=lter)
summary(errorCanopyFit, standardized=T)

error_pred_both <- as.data.frame(cbind(lavPredict(errorBothFit), inspect(errorCanopyFit, "data")))
names(error_pred_both) <- c("true_spring_biomass", "true_summer_biomass", "landsat_spring_biomass", "summer_kelp_biomass")

qplot(true_spring_biomass, true_summer_biomass, data=error_pred_both) +
  theme_bw(base_size=17)