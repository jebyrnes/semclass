#--------------------------------------------------------
# STRUCTURAL EQUATION MODELING: LOCAL ESTIMATION ###
#
# Author: Jon Lefcheck and Jarrett Byrnes
# Last updated: 2019-02-10
#--------------------------------------------------------

library(dagitty)

g <- dagitty("dag{
  pesticide -> caprellids <- macroalgae
  pesticide -> gammarids <- macroalgae
  caprellids -> epiphytes <- gammarids
  macroalgae -> epiphytes <- eelgrass
  eelgrass -> caprellids
  eelgrass -> gammarids
             }")

impliedConditionalIndependencies(g)


#adding exogenous covariances
g2 <- dagitty("dag{
  pesticide -> caprellids <- macroalgae
  pesticide -> gammarids <- macroalgae
  caprellids -> epiphytes <- gammarids
  macroalgae -> epiphytes <- eelgrass
  eelgrass -> caprellids
  eelgrass -> gammarids
  macroalgae <-> pesticide
  eelgrass <-> macroalgae
  eelgrass <-> pesticide
             }")
impliedConditionalIndependencies(g2)


#with lavaan

mod <- "
  epiphytes ~ eelgrass + macroalgae + caprellids + gammarids  
  caprellids ~ eelgrass + macroalgae + pesticide
  gammarids ~ eelgrass + macroalgae + pesticide
"

g3 <- lavaanToGraph(lavaan::lavaanify(mod))
impliedConditionalIndependencies(g3)


library(piecewiseSEM)

# Read in data
keeley <- read.csv("../data/Keeley_rawdata_select4.csv")

head(keeley)

### Fit Individual Relationships

abiotic_mod <- lm(abiotic ~ distance, data = keeley)

hetero_mod <- lm(hetero ~ distance, data = keeley)

rich_mod <- lm(rich ~ abiotic + hetero, data = keeley)


### Create list of structured equations
keeley.sem <- psem(

  abiotic_mod,

  hetero_mod,

  rich_mod,

  data = keeley

)

keeley.sem

### Evaluate fit

# Get the basis set
basisSet(keeley.sem)

# Conduct d-sep tests
claim1 <- lm(rich ~ distance + abiotic + hetero, keeley)

coefs(claim1)

claim2 <- lm(hetero ~ abiotic + distance, keeley)

coefs(claim2)

# Compute Fisher's C & compare to Chi-square distribution
(C <- -2 * (log(coefs(claim1)[1, 7]) + log(coefs(claim2)[1, 7])))

1 - pchisq(C, 2 * 2)

####

# Automagically conduct d-sep tests
keeley.dsep <- dSep(keeley.sem)

# By default the conditioning variables are hidden, but we can show them
dSep(keeley.sem, conditioning = TRUE)

# Get Fisher's C
fisherC(keeley.sem)

# Looks like one path (rich ~ distance) is HIGHLY significant!

# Add significant path back into model
keeley.sem2 <- update(keeley.sem, rich ~ abiotic + hetero + distance)

dSep(keeley.sem2)

fisherC(keeley.sem2)

# Get coefficients
coefs(keeley.sem2)

# Return intercepts as well
coefs(keeley.sem2, intercepts = T)

# Get R-squared
rsquared(keeley.sem2)

# Get all summary information
summary(keeley.sem2)

# Evaluate individual model assumptions & fits

# Graphical evaluation
par(mfrow = c(2, 2))

lapply(keeley.sem2[1:3], plot, which = 1)

### Exercise for fireseverity ####
mod_firesev <- lm(firesev ~ age, data=keeley)
mod_cover <- lm(cover ~ firesev, data=keeley)

firesev_fullmed_model <- psem(
  mod_firesev,
  mod_cover,
  keeley
)

dSep(firesev_fullmed_model)

coefs(firesev_fullmed_model)

rsquared(firesev_fullmed_model)

### Correlated errors

keeley.sem3 <- psem(

  lm(abiotic ~ distance, data = keeley),

  lm(hetero ~ distance, data = keeley),

  lm(rich ~ distance + hetero, data = keeley),

  rich %~~% abiotic

)

summary(keeley.sem3)

### Compare using AIC and LRT

keeley.sem4 <- psem(

  lm(hetero ~ distance, data = keeley),

  lm(rich ~ distance + hetero, data = keeley),
  
  abiotic ~ 1

)

anova(keeley.sem2, keeley.sem4)

AIC(keeley.sem2)

AIC(keeley.sem4)

### Get partial correlation plot
dev.off()

# Plot raw data
plot(rich ~ distance, keeley)

# Add trendline & 95% CIs
newdata.raw <- data.frame(distance =
                        seq(min(keeley[, 1], na.rm = TRUE) - min(keeley[, 1], na.rm = TRUE) * 0.1,
                            max(keeley[, 1], na.rm = TRUE) + min(keeley[, 1], na.rm = TRUE) * 0.1,
                            length.out = nrow(keeley) * 2) )


mod.raw <- lm(rich ~ distance, data = keeley)

pred.raw <- predict(mod.raw, newdata.raw, interval = "confidence", level = 0.95)

abline(mod.raw, col = "red", lwd = 2)

lines(newdata.raw[, 1], pred.raw[, 2], col = "red", lwd = 1.8, lty = 2)

lines(newdata.raw[, 1], pred.raw[, 3], col = "red", lwd = 1.8, lty = 2)

# Plot residuals
resids <- partialResid(rich ~ distance, keeley.sem2)


# Add trendline & 95% CIs
newdata <- data.frame(xresid =
                        seq(min(resids[, 2], na.rm = TRUE) - min(resids[, 2], na.rm = TRUE) * 0.1,
                            max(resids[, 2], na.rm = TRUE) + min(resids[, 2], na.rm = TRUE) * 0.1,
                            length.out = nrow(resids) * 2) )


#Partial Resid Plot
mod.resid <- lm(yresid ~ xresid, data = resids)

plot(yresid ~ xresid, xlab = "distance | others", ylab = "rich | others", data = resids)

abline(mod.resid, col="red", lwd=2)

#CIs
pred <- predict(mod.resid, newdata, interval = "confidence", level = 0.95)


lines(newdata[, 1], pred[, 2], col = "red", lwd = 1.8, lty = 2)

lines(newdata[, 1], pred[, 3], col = "red", lwd = 1.8, lty = 2)


#crPlots

plot(keeley$distance, pr$yresid, 
     xlab = "distance ", ylab = "rich | others")
cr.mod <- lm(yresid ~ distance, 
             data = data.frame(yresid=resids$yresid, distance=keeley$distance))
abline(cr.mod, col="red", lwd=2)

abline(mod.resid, col="red", lwd=2)


library(ggplot2)
ggplot(keeley, aes(x=keeley$distance, y = resids$yresid)) +
  geom_point() +
  stat_smooth(method="lm") +
  xlab("distance") + ylab("rich | others") +
  theme_bw(base_size=17)
