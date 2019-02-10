######################################################################
##### Code for Day 1 of SEM Workshop
##### covering ML fitting of SEMs of observed variables
#####
##### Jarrett E.K. Byrnes
#####
##### Last Modified 3/3/18
######################################################################


#Load your Data File
keeley<-read.csv("../data/Keeley_rawdata_select4.csv")

#load lavaan
library(lavaan)
library(lavaanPlot)


##############################
###Intro to Lavaan ####
##############################

#regression syntax versus lavaan
aLM<-lm(cover ~ age, data=keeley)

aSEM<-sem('cover ~ age', data=keeley)

#output from both
summary(aSEM)

summary(aLM)

lavaanPlot(model = aSEM, coefs = TRUE)
lavaanPlot(model = aSEM, coefs = TRUE,
           stand=TRUE)

#we want out intercept!
aMeanSEM<-sem('cover ~ age', data=keeley, meanstructure=T)

summary(aMeanSEM)

lavaanPlot(model = aMeanSEM, coefs = TRUE)


#standardized coefficients
standardizedSolution(aSEM)

summary(aSEM, standardized=T, rsq=T)

lavaanPlot(model = aMeanSEM, coefs = TRUE,
           stand = TRUE)


####################################
# SEM with multiple paths ####
####################################

partialMedModel<-' firesev ~ age
                cover ~ firesev + age'

partialMedSEM<-sem(partialMedModel, 
                   data=keeley)

#look at coefficients
summary(partialMedSEM, rsquare=T, standardized=T)


lavaanPlot(model = partialMedSEM, coefs = TRUE,
           stand = TRUE, 
           graph_options = list(layout = "circo"),
           sig = 0.05)


#################################
#Direct and Indirect effects ####
#################################

partialMedModelInd <-'

  #model
  firesev ~ af*age
  cover ~ fc*firesev + ac*age

  #Derived Calcuations
  indirect := af*fc
  total := ac + (af*fc)
'


partialMedSEMInd<-sem(partialMedModelInd, 
                   data=keeley)

summary(partialMedSEMInd)

standardizedSolution(partialMedSEMInd)

########################
#plot an sem ####
########################
library(lavaanPlot)

lavaanPlot(model = partialMedSEM, coefs = TRUE,
           stand = TRUE, 
           graph_options = list(layout = "circo"))


########################
#Exercise 1         ####
########################


#The Richness Partial Mediation Model
distModel <- 'rich ~ distance + abiotic + hetero
hetero ~ distance
abiotic ~ distance'

distFit <- sem(distModel, data=keeley)

summary(distFit, standardized=T, rsquare=T)
standardizedSolution(distFit)

distModelEff <- '
rich ~ dr*distance + ar*abiotic + hr*hetero
hetero ~ dh*distance
abiotic ~ da*distance

#The effects
direct := dr
indirect := dh*hr + da*ar
total := direct + indirect
'

distFitEff <- sem(distModelEff, data=keeley)

standardizedSolution(distFitEff)

########################
# Fixing Paths      ####
########################

#what if we know better
zeroMedModel<-'firesev ~ 0*age
              cover ~ 0*firesev + age'

zeroMedFit<-sem(zeroMedModel, data=keeley)
summary(zeroMedFit, rsquare=T)
standardizedSolution(zeroMedFit)


#Or use this intercept approach
zeroMedModel2<-'
firesev ~ 1
cover ~ age
'

zeroMedFit2<-sem(zeroMedModel2, data=keeley)
summary(zeroMedFit2, rsquare=T)
standardizedSolution(zeroMedFit2)

inspect(aSEM, "obs")
inspect(zeroMedFit, "obs")

######################
#Correlated Error ####
######################
corModel <-'firesev ~ age
            cover ~ age
            cover ~~ firesev'

corFit <- sem(corModel, data=keeley)
standardizedSolution(corFit)

####################
#Final Exercise ####
####################

#Part 1

corErrorModel <- '
  rich ~ distance + abiotic + hetero
  hetero ~ distance
  abiotic ~ distance

  abiotic ~~ hetero
'

corFit<-sem(corErrorModel, data=keeley)
summary(corFit, rsquare=T)
standardizedSolution(corFit)

#Part 2
oneDistModel <- 'rich ~ 1*distance + abiotic + hetero
              hetero ~ distance
              abiotic ~ distance'
oneFit<-sem(oneDistModel, data=keeley)
summary(oneFit, rsquare=T)
standardizedSolution(oneFit)
