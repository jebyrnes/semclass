#'-----------------------------------------------------------------------
#' STRUCTURAL EQUATION MODELING: Nonlinear Nonnormal Estimation ###
#'
#' Author: Jon Lefcheck & Jarrett Byrnes
#' Last updated: 3 March 2017
#'-----------------------------------------------------------------------

library(piecewiseSEM)
library(lavaan)


## Non-linearities (Cardinale 2009)

# Load data
cardinale <- read.csv("../data/cardinale.csv")

head(cardinale)

# Take log of N and N^2
cardinale$logN <- log10(cardinale$N + 1e-6)

cardinale$logN2 <- cardinale$logN ^ 2

# Take log of chl (standing biomass)
cardinale$logChl <- log10(cardinale$Chl)

# Create SEM
cardinale.sem <- psem(
  lm(SA ~ logN + logN2 + SR, data = cardinale),
  lm(logChl ~ SA + logN + logN2, data = cardinale),
  cardinale
  )

fisherC(cardinale.sem)

# Evaluate fit
summary(cardinale.sem)
coefs(cardinale.sem)

with(cardinale, cor(logN, logN2))


# Could model fits be influenced by high correlation between squared and unsquared N terms?

# Center polynomial to reduce collinearity
cardinale$logN.cen = scale(cardinale$logN, scale = F)

cardinale$logN2.cen = scale(cardinale$logN, scale = F) ^ 2

cor(cardinale$logN.cen, cardinale$logN2.cen)

# Re-fit SEM using centered predictors
cardinale.sem2 <- psem(
  lm(SA ~ logN.cen + logN2.cen + SR, data = cardinale),
  lm(logChl ~ SA + logN.cen + logN2.cen, data = cardinale),
  cardinale
)


# Evaluate fit
coefs(cardinale.sem2)
summary(cardinale.sem2)

dSep(cardinale.sem2)

# with lavaan

cardinale_mod <-'
SA ~ logN.cen + logN2.cen + SR
logChl ~ SA + logN.cen+ logN.cen + logN2.cen
logN.cen ~~ logN2.cen
'
cardinale_fit <- sem(cardinale_mod, data=cardinale, fixed.x=F)

summary(cardinale_fit, std=TRUE)

# Read in data
keeley <- read.csv("../data/Keeley_rawdata_select4.csv")


### Create list of structured equations
keeley_int <- psem(

  lm(cover ~ age*firesev, data = keeley),

  lm(firesev ~ age, data = keeley),


  data = keeley

)

fisherC(keeley_int)
coefs(keeley_int)

with(keeley, cor(age, age*firesev))

## with centering

keeley$fire_cent <- scale(keeley$firesev, scale=FALSE)
keeley$age_cent <- scale(keeley$age, scale=FALSE)

  
  
keeley_int_cent <- psem(
  
  lm(cover ~ age_cent *fire_cent, data = keeley),
  
  lm(fire_cent ~ age_cent, data = keeley),
  
  
  data = keeley
  
)

coefs(keeley_int_cent)

#and with lavaan
partialMedModel_int<-' 
  firesev ~ age

  cover ~ firesev + age + firesev:age
'

partialMedFit_int <- sem(partialMedModel_int, data=keeley)

summary(partialMedFit_int)
modificationIndices(partialMedFit_int)


partialMedModel_int_scale<-' 
  fire_cent ~ age_cent

  cover ~ fire_cent + age_cent + fire_cent:age_cent
'

partialMedFit_int_scale <- sem(partialMedModel_int_scale, data=keeley)


summary(partialMedFit_int_scale)
### Fit Anderson SEM with GLM ####
library(ggplot2)

# Load data
anderson <- read.csv("../data/anderson.csv")

qplot(biomass.kg, hotspotYN, data=anderson) +
  theme_bw(base_size=17) +
  stat_smooth(method="glm", method.args=list(family=binomial), color="red", lwd=2)
  

# Generate SEM
anderson.sem <- psem(
  
  lm(leafN ~ biomass.kg, anderson),
  
  glm(hotspotYN ~ leafN + biomass.kg + landscape, family = "binomial", anderson),
  
  data = anderson
  
)

# Recover summary statistics and coefficients
summary(anderson.sem)
rsquared(anderson.sem)

### A WARNING ####

# Create fake data
set.seed(1)

data <- data.frame(
  x = runif(100),
  y1 = runif(100),
  y2 = rpois(100, 1),
  y3 = runif(100)
)

# Create SEM with GLM
modelList <- psem(
  lm(y1 ~ x, data),
  glm(y2 ~ x, "poisson", data),
  lm(y3 ~ y1 + y2, data),
  data
)

# Show that y2 ~ y1 is the same as y2 ~ y1 for LM
mody1.y2 <- lm(y1 ~ y2 + x, data)

mody2.y1 <- lm(y2 ~ y1 + x, data)

summary(mody1.y2)$coefficients[2, 4]

summary(mody2.y1)$coefficients[2, 4]

# Show that y2 ~ y1 is not the same as y2 ~ y1 for GLM
mody1.y2 <- lm(y1 ~ y2 + x, data)

mody2.y1.glm <- glm(y2 ~ y1 + x, "poisson", data)

summary(mody1.y2)$coefficients[2, 4]

summary(mody2.y1.glm)$coefficients[2, 4]

# Run summary
summary(modelList)

# Address conflict using conserve = T
summary(modelList, conserve = T)

dSep(modelList, conserve = T)

summary(mody1.y2)$coefficients[2, 4]

summary(mody2.y1.glm)$coefficients[2, 4]

# Address conflict using direction = c()
summary(modelList, direction = c("y2 <- y1"))

dSep(modelList, direction = c("y2 <- y1"))

dSep(modelList, direction = c("y1 <- y2"))

# Address conflict using correlated errors
modelList2 <- update(modelList, y2 %~~% y1)

dSep(modelList2)

# Now with the anderson data
fisherC(anderson.sem, conserve=TRUE)

## GLMs and standardized coefficients ####

# Subset just GLM
anderson.glm <- anderson.sem[[2]]

# Get vector of Beta values
Betas <- coefs(anderson.sem)[2:4, 3]

# Compute predicted values on the link scale
preds <- predict(anderson.glm, type = "link")
preds_response <- predict(anderson.glm, type = "response")

# Compute sd of predicted  values using theoretical error variance
sd.y.LT <- sqrt(var(preds) + pi^2/3)

# Compute sd of predictors
sd.x <- sapply(anderson[, c("leafN", "biomass.kg", "landscape")], sd)

# Scale Betas by sd.x / sd.y
(Betas.LT <- Betas * sd.x / sd.y.LT)

# Compute sd of predicted values using observed variance
R <- cor(anderson$hotspotYN, predict(anderson.glm, type = "response"))

sd.y.OE <- sqrt(var(response)) / R

# Get standardized coefficient
(Betas.OE <- Betas * sd.x / sd.y.OE)

# Look at partial correlations
cor(anderson[, c("hotspotYN", "leafN", "biomass.kg", "landscape")])

# Get standardized effect of biomass.kg on leafN
Beta.leafN <- coefs(anderson.sem)$Std.Estimate[1]

# Compute indirect effect based on LT then compare to direct effect
Beta.leafN * Betas.LT[1]; Betas.LT[2]

# Compute indirect effect based on LT then compare to direct effect
Beta.leafN * Betas.OE[1]; Betas.OE[2]

# LOGITS and std coefs
logitcurve <- data.frame(y = rbinom(100, prob=seq(0,1,length.out=100), size=4)/4, x=1:100)

logitmod <- glm(y ~ x, data=logitcurve,
                weights=rep(4,nrow(logitcurve)),
                family="binomial")
predfun <- function(x) predict(logitmod, newdata=data.frame(x=x), type="response")

#Plot
par(mfrow=c(1,2))
plot(y ~ x, data=logitcurve, pch=19,
     xlab="Body Size", ylab="Probability of Survival")
curve(predfun, lwd=3, col="red", from=1, to=100, add=TRUE)

plot(boot::inv.logit(y) ~ x, data=logitcurve, pch=19,
       xlab="Body Size", 
     ylab="Inv. Logit Probability of Survival")
abline(lm(boot::inv.logit(y) ~ x, data=logitcurve),
       lwd=3, col="red")

par(mfrow=c(1,1))

#anderson std coefs
#The default is Menard.OE
coefs(anderson.sem, standardize.type = "Menard.OE")
coefs(anderson.sem, standardize.type = "latent.linear")

