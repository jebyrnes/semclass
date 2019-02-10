### STRUCTURAL EQUATION MODELING: AUTOCORRELATION ###

### Author: Jon Lefcheck
### Last updated: 24 October 2017

#install.packages("ape")
library(ape)
library(car)
library(lavaan)
library(MASS)
library(nlme)
library(piecewiseSEM)

#install.packages("spdep")
library(spdep)

#install.packages("ncf")
library(ncf)

#install.packages("sp")
library(sp)


### SPATIAL WEIGHTS EXAMPLE

# Read in data
boreal <- read.csv("../Data/boreal.csv")

head(boreal)
#Plot it
library(ggplot2)
qplot(x, y, data=boreal, size=wet, color=NDVI) +
  theme_bw(base_size=18) + 
  scale_size_continuous("Index of Wetness", range=c(0,10)) + 
  scale_color_gradient("NDVI", low="lightgreen", high="darkgreen")

qplot(x, y, data=boreal, size=richness, color=NDVI) +
  theme_bw() + 
  scale_size_continuous("Species Richness", range=c(0,10)) + 
  scale_color_gradient("NDVI", low="lightgreen", high="darkgreen")

qplot(x, y, data=boreal, size=temp, color=NDVI) +
  theme_bw() + 
  scale_size_continuous("Days Freezing", range=c(0,10)) + 
  scale_color_gradient("NDVI", low="lightgreen", high="darkgreen")


# Run regression models (no spatial weights)

rich_lm <- lm(richness ~ temp, data = boreal)

ndvi_lm <- lm(NDVI ~ richness + temp + wet, data=boreal)

boreal.sem <- psem(

  rich_lm,

  ndvi_lm,

  temp %~~% wet,

  data = boreal

)

summary(boreal.sem)

# Extract residuals
source("./residuals.R")

res <- residuals.psem(boreal.sem)
boreal <- cbind(boreal, res)

#visualize

qplot(x, y, data=boreal, color=NDVI_residuals, size=I(5)) +
  theme_bw(base_size=17) + 
  scale_color_gradient("NDVI Residual", low="blue", high="yellow")

qplot(x, y, data=boreal, size=NDVI_residuals, color=factor(NDVI_residuals>0)) +
  theme_bw() + scale_size_continuous("Absolute Value of Residual", range=c(0,5)) + scale_color_discrete("Residual > 0?")

qplot(x, y, data=boreal, color=richness_residuals, size=I(5)) +
  theme_bw(base_size=17) + 
  scale_color_gradient("Richness Residual", low="blue", high="yellow")

#Moran's I
library(ape)
distMat <- as.matrix(dist(cbind(boreal$x, boreal$y)))

distsInv <- 1/distMat
diag(distsInv) <- 0

mi.ndvi <- Moran.I(boreal$NDVI_residuals, distsInv)
mi.ndvi


Moran.I(boreal$richness_residuals, distsInv)

#Visualize variogram
library(nlme)
ndvi_gls<- gls(NDVI ~ richness + temp + wet, data=boreal)

plot(Variogram(ndvi_gls, form=~x+y, 
               robust=T, maxDist=2000, 
               resType="normalized"))

#Fit using spatial autocorrelation
spaceCor <- corExp(form =~ x+y, nugget=T)

ndvi_gls_space<- gls(NDVI ~ richness + temp + wet, 
                     correlation = spaceCor,
                     data=boreal)

rich_gls_space <- gls(richness ~ temp,
                     correlation = spaceCor,
                     data = boreal)

boreal_space.sem <- psem(
  
  ndvi_gls_space,
  
  rich_gls_space,
  
  temp %~~% wet,
  
  data = boreal
  
)

#evaluate
dSep(boreal_space.sem)

coefs(boreal_space.sem)

### Nearest Neighbors approach

# Determine nearest neighbors
boreal_sp <- boreal
coordinates(boreal_sp) <- ~x+y

nb <- tri2nb(boreal_sp)

plot(nb, coordinates(boreal_sp))

# Run regression models (with spatial weights)
spatial_weights <- nb2listw(nb)

rich_lag <- lagsarlm(richness ~ temp, 
                     data = boreal_sp,
                     listw = spatial_weights, 
                     tol.solve = 1e-11)

ndvi_lag <- lagsarlm(NDVI ~ richness + temp + wet, 
                     data = boreal_sp, 
                     listw = spatial_weights, 
                     tol.solve = 1e-11)

boreal_space_lag.sem <- psem(

  rich_lag,

  ndvi_lag,

  temp %~~% wet,

  data = boreal_sp

)

summary(boreal.sem2)

## Example
squid <- read.csv("../Data/goodsquid.csv")

head(squid)

qplot(SqYear, Squid, data = squid)

ggplot(squid, aes(LonW, LatN, color = log(SquidCount+1))) + 
  geom_point() +
  scale_color_gradient(low="yellow", high="purple")

ggplot(squid %>% filter(SquidCount>0), aes(LonW, LatN, color = log(SquidCount+1))) + 
  geom_point() +
  scale_color_gradient(low="yellow", high="purple")

hake_mod <- gls(Hake ~ UIwin + DepthOMZ + Distkm, data = squid)

plot(squid$SqYear, residuals(hake_mod))


plot(Variogram(hake_mod, form=~ LatN + LonW, 
               robust=T, 
               resType="normalized", nint=500))

corSquid <- corExp( form = ~ LatN + LonW | SqYear/SqMonth, nugget = T)

hake_mod_space <- gls(Hake ~ UIwin + DepthOMZ, data = squid,
                      correlation = corSquid)

summary(hake_mod_space)


plot(Variogram(hake_mod_space, form=~ LatN + LonW, 
               robust=T, maxDist=1, 
               resType="normalized", nint=500))
##### Data processing
# 
# squid <- read.csv("../Data/squid_stewart.csv")
#load("../../range-expansion-modeling-master/data/mysquid.RData")
#write.csv(mysquid, "../data/stewart_squid.csv", row.names=FALSE)
# #SqMonth: calculated 12-month year defined from March (3) - Feb (14) (MBARI)
# squid$month <- squid$SqMonth
# squid$month[which(squid$month==13)] <- 1
# squid$month[which(squid$month==14)] <- 2
# 
# squid$date = with(squid, 
#                   lubridate::parse_date_time(paste(month, SqYear, sep=","),
#                                              orders="my"))
# squid %>% group_by(date, LatN, LonW) %>% count() %>% ungroup() %>% filter(n>1) -> badSquid
# # badSquid
# # # A tibble: 591 x 4
# # date   LatN     LonW     n
# # <dttm>  <dbl>    <dbl> <int>
# #   1 1998-04-30 36.699 -122.031     2
# # 2 1999-01-31 36.583 -122.516     2
# # 3 2002-04-30 36.330 -122.900     2
# # 4 2004-05-31 36.330 -122.900     2
# # 5 2005-08-31 36.330 -122.899     2
# # 6 2005-08-31 36.736 -122.047     2
# # 7 2006-03-31 36.330 -122.900     2
# # 8 2006-04-30 36.712 -122.189     2
# # 9 2008-05-31 36.802 -121.994     2
# 
# goodSquid <- squid %>% group_by(date, LatN, LonW) %>%
#   slice(1L) %>%
#   ungroup()
# head(goodSquid)
# nrow(goodSquid)
# nrow(squid)
# write.csv(goodSquid, file="../data/goodsquid.csv", row.names=FALSE)
