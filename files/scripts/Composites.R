######################################################################
##### Code for Day 4 of SEM Workshop
##### covering Composite Variables 
#####
##### Jarrett E.K. Byrnes
#####
##### Last Modified 10/21/17
######################################################################


#####################################################################################
# Composite Variables for Nonlinearities
###################################################################################################


#Load your Data File
keeley<-read.csv("../data/Keeley_rawdata_select4.csv")


#is this a nonlinear relationship?
plot(rich ~ cover, data=keeley, cex.lab=1.5, cex.axis=1.3)

#compare models
linear<-lm(rich ~ cover, data=keeley)
nonlinear<-lm(rich ~ cover+I(cover^2), data=keeley)

library(AICcmodavg)
aictab(list(linear, nonlinear), c("linear", "squared"))

#plot the nonlinear result
curve(coef(nonlinear)[1] + coef(nonlinear)[2]*x + coef(nonlinear)[3]*x^2, add=T, col="red", lwd=4)
curve(coef(linear)[1] + coef(linear)[2]*x, add=T, col="blue", lwd=4)

#Step 1: Create a new nonlinear variable in the data
keeley<-within(keeley, coverSQ<-cover^2)

#Now, we can fit this with no composite
noCompModel <- 'rich ~ cover + coverSQ'

noCompFit <- sem(noCompModel, data=keeley)
summary(noCompFit, rsquare=T)

#Now, add a composite

compModel<-'

					coverEffect <~ cover + coverSQ

				  rich ~ 1*coverEffect'

compFit <- sem(compModel, data=keeley)

summary(compFit)

###
#Exercise: Abiotic
###

#how about abotic and fire severity?
plot(firesev ~ abiotic, data=keeley)

fireLinear<-lm(firesev ~ abiotic, data=keeley)
fireNonlinear<-lm(firesev ~ abiotic + I(abiotic^2), data=keeley)

aictab(list(fireLinear, fireNonlinear), c("linear", "squared"))

curve(coef(fireNonlinear)[1] + coef(fireNonlinear)[2]*x + coef(fireNonlinear)[3]*x^2, add=T, col="red", lwd=4)

keeley$abioticSQ <- keeley$abiotic^2



abioticCompModelBad<-'
                   abioticEffect <~ 0.4 * abiotic + abioticSQ

                  firesev ~ abioticEffect'

abioticCompModel<-'
                   abioticEffect <~ abiotic + abioticSQ

                  firesev ~ 1*abioticEffect'

abioticCompFit <- sem(abioticCompModel, data=keeley)

summary(abioticCompFit)

####
# Endogenous composite
####


endoCompModel<-'
                coverEffect <~ cover + coverSQ

                cover ~~ coverSQ
                age ~~ coverSQ

                cover ~ age
                rich ~ age + 1*coverEffect'


endoCompFit <- sem(endoCompModel, data=keeley, fixed.x=F)

summary(endoCompFit)



###########
# Composite in a piecewise approach
###########



cardModel<-'
SA ~ logN + logN2Cen + SR
logChl ~ SA + logN
'

cardModelComposite <- '
N_effect <~ logN + logN2Cen
SA ~ 1*N_effect + SR
logChl ~ SA + logN
'
carCompositeFit <- sem(cardModelComposite, data=cards)
standardizedSolution(carCompositeFit)

#with lavaan

fireCompModel<-'
  #composite definitions
	coverEffect <~ cover + coverSQ 

  #Covariance bookkeeping
  cover ~~ coverSQ
  coverSQ ~~ firesev

  #model
  cover ~ firesev
  rich ~ 1*coverEffect + firesev

'
fireCompFit <- sem(fireCompModel, data=keeley)


#First, fit the observed only model on the composite piece
rich_obs_mod <- lm(rich ~ cover + coverSQ + firesev, data=keeley)

#Now extract a composite
keeley$cover_eff <- with(keeley, coef(rich_obs_mod)[2]*cover + coef(rich_obs_mod)[3]*cover^2)

#Second, make a loadings relationship
comp_loadings_mod <- lm(cover_eff ~ cover + coverSQ, data=keeley)

#Third, refit the model with the composite
rich_comp_mod <- lm(rich ~ offset(1*cover_eff) + offset(-2.161*firesev), data=keeley)

#Now, put it all together
cover_mod <- lm(cover ~ firesev, data=keeley)

#Roll it into a pSEM
fire_comp_psem <- psem(
  comp_loadings_mod,
  cover_mod,
  rich_comp_mod,
  coverSQ %~~% firesev,
  coverSQ %~~% cover,
  data = keeley
)


coefs(fire_comp_psem)
basisSet(fire_comp_psem)
dSep(fire_comp_psem)
fisherC(fire_comp_psem)
summary(fire_comp_psem)
