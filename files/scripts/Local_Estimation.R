#--------------------------------------------------------
# STRUCTURAL EQUATION MODELING: LOCAL ESTIMATION ###
#
# Author: Jon Lefcheck and Jarrett Byrnes
# Last updated: 2019-02-10
#--------------------------------------------------------

library(dagitty)

g <- dagitty("dag{
  x1 -> y1
  x1 -> y2
  x1 -> y3
  y3 -> y4
  y2 -> y4
  y1 -> y4
  y2 -> y1
             }")

#conditional independence relationships
impliedConditionalIndependencies(g)



# My Model
forest_mod <-  dagitty("dag{
  	waves -> kelp -> algae  
	  algae -> inverts
  	waves -> algae
}")


plot(graphLayout(forest_mod))

adjustmentSets(forest_mod, 
               exposure = "kelp", 
               outcome = "inverts")

#conditional independence relationships
impliedConditionalIndependencies(forest_mod)


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

#plot the model
plot(keeley.sem2)

#make a better plot
keeley.plot2 <- plot(keeley.sem2, return=TRUE)

#what's in there?
get_node_df(keeley.plot2)

keeley.plot2 %>%
  set_node_attrs(node_attr = x, values = c(2.5,2.5,4,1)) %>%
  set_node_attrs(node_attr = y, values = c(3,1,2,2)) %>%
  render_graph()

#or...
plot(keeley.sem2,
     node_attrs = list(x = c(2.5,2.5,4,1),
                       y = c(3,1,2,2),
                       shape = "rectangle", fillcolor = "white"))

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

# Plot raw data
plot(rich ~ distance, keeley)

#psem objects are lists of models
keeley.sem2[[3]]

#Added Variable Plots
library(car)

avPlot(keeley.sem2[[3]], variable = "distance")


#CR Plots
crPlot(keeley.sem2[[3]], variable = "distance",
       smooth = FALSE)

#visreg
library(visreg)

visreg(keeley.sem2[[3]], 
       xvar = "distance")

vr <- visreg(keeley.sem2[[3]], 
             xvar = "distance")

head(vr$res)


visreg(keeley.sem2[[3]], 
       xvar = "distance",
       gg=TRUE) + ggthemes::theme_par()

#Final Exercise
fullMed <- psem(
  lm(cover ~ firesev, data = keeley),
  lm(firesev ~ age, data = keeley),
  keeley)
)

noMed <- psem(
  lm(firesev ~ age, data = keeley),
  lm(cover ~ age, data = keeley),
  data = keeley)

anova(fullMed, noMed)

#visualize
partialMed <- psem(
  lm(firesev ~ age, data = keeley),
  lm(cover ~ age + firesev, data = keeley),
  data = keeley)
  
visreg(partialMed[[2]], xvar = "firesev")
