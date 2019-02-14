### STRUCTURAL EQUATION MODELING: MIXED MODELS ###

### Author: Jon Lefcheck & Jarrett Byrnes
### Last updated: 30 October 2017

library(ggplot2)
library(nlme)
library(lme4)
library(lavaan)
library(piecewiseSEM)

### COMPARISON OF MIXED MODEL STRUCTURES

## Varying slope, varying intercept models

# Create fake levels
set.seed(123)

data <- data.frame(
  level = rep(LETTERS[1:3], each = 50),
  x = runif(150, 0, 5)
)

# Generate different sample relationships within each level
data <- cbind(data,

             do.call(rbind, lapply(1:length(unique(data$level)), function(i) {

               datf <- subset(data, level == unique(level)[i])

               datf$y <- (i^2 + (2 *i) * i^2) + (sqrt(i) + i^2) * datf$x + runif(50, 0, 50)

               return(datf[, "y", drop = F])

             } ) )

)

## Fit varying intercept, fixed slope
vint.mod <- lme(y ~ x, random = ~ 1 | level, data)

# Look at coefficients for each level of the random factor
coef(vint.mod)

# Generate new data
newdata <- expand.grid(
  level = unique(data$level),
  x = seq(0, max(data$x), 0.01)
)

# Store predictions
newdata$vint.pred = predict(vint.mod, newdata)

# Plot results
ggplot() +
  geom_point(data = data, aes(x = x, y = y), col = "grey10") +
  geom_line(data = newdata, aes(x = x), y = predict(vint.mod, newdata, level = 0), col = "black", lwd = 1) +
  geom_line(data = newdata, aes(x = x, y = vint.pred, group = level, col = level), lwd = 1) +
  theme_bw(base_size = 18)

## Fit varying intercept, varying slope

vint.vslope.mod <- lme(y ~ x, random = ~ x | level, control = lmeControl(opt = "optim"), data)

# Look at coefficients for each level of the random factor
coef(vint.vslope.mod)

# Store predictions
newdata$vint.vslope.pred <- predict(vint.vslope.mod, newdata)

# Plot results
ggplot() +
  geom_point(data = data, aes(x = x, y = y), col = "grey10") +
  geom_line(data = newdata, aes(x = x), y = predict(vint.vslope.mod, newdata, level = 0), col = "black", lwd = 1) +
  geom_line(data = newdata, aes(x = x, y = vint.vslope.pred, group = level, col = level), lwd = 1) +
  theme_bw(base_size = 18)

## Fit fixed intercept, varying slope

vslope.mod <- lme(y ~ x, random = ~ 0 + x | level, data)

# Look at coefficients for each level of the random factor
coef(vslope.mod)

# Store predictions
newdata$vslope.pred <- predict(vslope.mod, newdata)

# Plot results
ggplot() +
  geom_point(data = data, aes(x = x, y = y), col = "grey10") +
  geom_line(data = newdata, aes(x = x), y = predict(vslope.mod, newdata, level = 0), col = "black", lwd = 1) +
  geom_line(data = newdata, aes(x = x, y = vslope.pred, group = level, col = level), lwd = 1) +
  theme_bw(base_size = 18)

## Nested models

names(data)[1] <- "level1"

# Add site variable
data$level2 <- rep(letters[1:2], each = 25)

# Create new model
nested.mod <- lme(y ~ x, random = ~ 1 | level1 / level2, control = lmeControl(opt = "optim"), data)

coef(nested.mod)










# Generate new data
newdata2 <- expand.grid(
  level1 = unique(data$level1),
  level2 = unique(data$level2),
  x = seq(0, max(data$x), 0.01)
)

# Store predictions
newdata2$nested.pred <- predict(nested.mod, newdata2, level = 0:1)

# Plot results

### R-SQUARED

## Comparison of OLS vs ML for linear regression

# Generate fake data
data <- data.frame(y = runif(100))

data$x <- data$y + runif(100, 0, 0.5)

# OLS
mod <- lm(y ~ x, data)

SSE <- sum((data$y - fitted(mod))^2)

SSTot <- sum((data$y - mean(data$y))^2)

(SS.R2 <- 1 - SSE / SSTot); summary(mod)$r.squared

# Likelihood
nullmod <- lm(y ~ 1, data)

L1 <- as.numeric(logLik(mod))

L0 <- as.numeric(logLik(nullmod))

L.R2 <- 1 - exp(-2/nrow(data) * (L1 - L0))

SS.R2; L.R2; summary(mod)$r.squared

# Refit for GLM
data$y <- rbinom(100, 1, 0.5)

mod.glm <- glm(y ~ x, "binomial", data)

SSE <- sum((data$y - fitted(mod.glm))^2)

SSTot <- sum((data$y - mean(data$y))^2)

SS.R2 <- 1 - SSE / SSTot

nullmod.glm <- glm(y ~ 1, "binomial", data)

L1 <- as.numeric(logLik(mod.glm))

L0 <- as.numeric(logLik(nullmod.glm))

L.R2 <- 1 - exp(-2/nrow(data) * (L1 - L0))

SS.R2; L.R2 # not identical

#R2 for mixed model
fm1 <- lme(distance ~ age, data = Orthodont) # random is ~ age

rsquared(fm1)

### MIXED MODEL SEM

# Load data
shipley <- read.csv("../Data/shipley.csv")


# dd_mod <- lme(DD ~ lat, random = ~1|site/tree, na.action = na.omit,
#     data = shipley)
# 
# date_mod <- lme(Date ~ DD, random = ~1|site/tree, na.action = na.omit,
#     data = shipley)
# 
# growth_mod <- lme(Growth ~ Date, random = ~1|site/tree, na.action = na.omit,
#     data = shipley)
# 
# live_mod <- glmer(Live ~ Growth + (1|site) + (1|tree),
#       family=binomial(link = "logit"), data = shipley) 

#missing values?
library(visdat)
vis_dat(shipley)

shipley <- na.omit(shipley)
vis_dat(shipley)


dd_mod <- lmer(DD ~ lat + (1|site/tree), 
               data = shipley)

date_mod <- lmer(Date ~ DD + (1|site/tree), 
                 data = shipley)

growth_mod <- lmer(Growth ~ Date + (1|site/tree), 
                   data = shipley)

live_mod <- glmer(Live ~ Growth + (1|site/tree),
                  family=binomial(link = "logit"), 
                  data = shipley)


# Create list of structural equations
shipley.sem <- psem(
  dd_mod,
  date_mod,
  growth_mod,
  live_mod,
  data = shipley

)

# Get summary
summary(shipley.sem, conserve = TRUE)

# Plot residuals vs fitted values
par(mfrow = c(2, 2))

lapply(shipley.sem, plot)

#residuals for a glmm
library(DHARMa)
sims <- simulateResiduals(shipley.sem[[4]])
plot(sims)

# Get indirect effect of latitude on survival
prod(coefs(shipley.sem)$Std.Estimate)

### Alternate model

# Create list of structural equations
shipley.sem2 <- psem(

  lme(DD ~ lat, random = ~1|site/tree, na.action = na.omit,
      data = shipley),

  lme(Date ~ DD, random = ~1|site/tree, na.action = na.omit,
      data = shipley),

  lme(Growth ~ DD, random = ~1|site/tree, na.action = na.omit,
      data = shipley),

  glmer(Live ~ Growth + (1|site) + (1|tree),
        family=binomial(link = "logit"), data = shipley),
  
  shipley

)

# Compare using AIC
AIC(shipley.sem2, shipley.sem) # definitely model 1!

# Re-fit using lavaan
shipley.model <- '
DD ~ lat
Date ~ DD
Growth ~ DD
Live~ Growth
'

shipley.lavaan <- sem(shipley.model, shipley)

summary(shipley.lavaan)

### Phytoplankton example

# Load data
durocher <- read.csv("../data/durocher.csv")

# names(durocher) <- c("Pond.ID", "Month", "Treatment", "Temp", "GPP", "Metabolic", "Biomass", "Zbio", "Richness")

head(durocher)

# Create model list
durocher.sem <- psem(

  lme(CR ~ Std.Temp, random = ~ 1 | Pond.ID, na.action = na.omit, durocher),

  lme(Prich ~ Std.Temp, random = ~ 1 | Pond.ID, na.action = na.omit, durocher),

  lme(Pbio ~ Prich, random = ~ 1 | Pond.ID, na.action = na.omit, durocher),

  lme(GPP ~ Pbio, random = ~ 1 | Pond.ID, na.action = na.omit, durocher)

)

# Evaluate fit
summary(durocher.sem)

# Remove NAs
durocher2 = durocher[complete.cases(durocher), ]

durocher.sem2 <- update(durocher.sem, data = durocher2)

summary(durocher.sem2)

# Play with random structure
durocher.sem3 <- psem(

  lme(CR ~ Std.Temp, random = ~ 1 | Pond.ID/Month, na.action = na.omit, durocher2),

  lme(Prich ~ Std.Temp, random = ~ 1 | Pond.ID/Month, na.action = na.omit, durocher2),

  lme(Pbio ~ Prich, random = ~ 1 | Pond.ID/Month, na.action = na.omit, durocher2),

  lme(GPP ~ Pbio, random = ~ 1 | Pond.ID/Month, na.action = na.omit, durocher2)

)

summary(durocher.sem3)

# By treatment
durocher.A.sem <- psem(

  lme(CR ~ Std.Temp, random = ~ 1 | Pond.ID/Month, na.action = na.omit, subset(durocher2, Treatment == "A")),

  lme(Prich ~ Std.Temp, random = ~ 1 | Pond.ID/Month, na.action = na.omit, subset(durocher2, Treatment == "A")),

  lme(Pbio ~ Prich, random = ~ 1 | Pond.ID/Month, na.action = na.omit, subset(durocher2, Treatment == "A")),

  lme(GPP ~ Pbio, random = ~ 1 | Pond.ID/Month, na.action = na.omit, subset(durocher2, Treatment == "A")),

  data = subset(durocher2, Treatment == "A")

)

summary(durocher.A.sem)

durocher.H.sem <- psem(

  lme(CR ~ Std.Temp, random = ~ 1 | Pond.ID/Month, na.action = na.omit, subset(durocher2, Treatment == "H")),

  lme(Prich ~ Std.Temp, random = ~ 1 | Pond.ID/Month, na.action = na.omit, subset(durocher2, Treatment == "H")),

  lme(Pbio ~ Prich, random = ~ 1 | Pond.ID/Month, na.action = na.omit, subset(durocher2, Treatment == "H")),

  lme(GPP ~ Pbio, random = ~ 1 | Pond.ID/Month, na.action = na.omit, subset(durocher2, Treatment == "H")),

  data = subset(durocher2, Treatment == "H")

)

summary(durocher.H.sem)

# How does one reproduce their analysis?!

# cards

## Non-linearities (Cardinale 2009)

# Load data
cardinale <- read.csv("../data/cardinale.csv")

head(cardinale)

# Take log of N and N^2
cardinale$logN <- log10(cardinale$N + 1e-6)

cardinale$logN2 <- cardinale$logN ^ 2

# Take log of chl (standing biomass)
cardinale$logChl <- log10(cardinale$Chl)

# Center polynomial to reduce collinearity
cardinale$logN.cen = scale(cardinale$logN, scale = F)

cardinale$logN2.cen = scale(cardinale$logN, scale = F) ^ 2


# Re-fit SEM using centered predictors
cardinale.sem2 <- psem(
  lm(SA ~ logN.cen + logN2.cen + SR, data = cardinale),
  lm(logChl ~ SA + logN.cen + logN2.cen, data = cardinale),
  cardinale
)


# Re-fit SEM using centered predictors
cardinale.mixed <- psem(
  lme(SA ~ logN.cen + logN2.cen + SR, 
      random = ~1|Stream, data = cardinale),
  
  lme(logChl ~ SA + logN.cen + logN2.cen, 
     random = ~1|Stream,  data = cardinale),
  
  data = cardinale
)

coefs(cardinale.mixed)

rsquared(cardinale.sem2)
rsquared(cardinale.mixed)

## Fully hierarchical Model
library(ggplot2)
theme_set(theme_bw(base_size=17))
ostra_site <- read.csv("../data/ostracod_sitelevel.csv")
ostra_plot <- read.csv("../data/ostracod_plotlevel.csv")


qplot(winter_temp, winter_phytoplankton, color=site, data = lowerlevel)


qplot(algal_cover, ostracod_abund, color=site, data = lowerlevel)
qplot(winter_phytoplankton, ostracod_abund, color=site, data = lowerlevel)
qplot(winter_temp, ostracod_abund, color=site, data = lowerlevel)

#site level
ostra_site_model <- psem(
  lm(winter_phytoplankton ~ winter_temp, data = ostra_site),
  
  data = ostra_site
)

#plot level
ostra_plot_model <- psem(
  lme(ostracod_abund ~ algal_cover + winter_phytoplankton,
      random = ~1|site, data = ostra_plot),
  
  data = ostra_plot
)

fisherC(ostra_site_model)
fisherC(ostra_plot_model)


basis_mod <-  lme(ostracod_abund ~ algal_cover + winter_phytoplankton + winter_temp,
                  random = ~1|site, data = ostra_plot)

summary(basis_mod)

fish_c <- 0 + 0 + -2*log(0.6275)

1 - pchisq(fish_c, df = 1)
