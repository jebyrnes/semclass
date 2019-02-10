######################
# Temporal Autocorrelation in SEMs 
#
# Jarrett Byrnes
#
# Last Updated 10/30/17
######################


library(car)
library(MASS)
library(nlme)
library(piecewiseSEM)

### Panel Model
spurge <- read.csv("../data/spurge00-01_data_logged.csv")

### TEMPORAL AUTOCORRELATION EXAMPLE
spurge_mod <- psem(
  lm(LNIG00 ~ LSPAD00, data=spurge),
  lm(LLAC00 ~ LSPAD00, data=spurge),
  lm(LNIG01 ~ LNIG00 + LSPAD00 + LLAC00, data = spurge),
  lm(LLAC01 ~ LLAC00 +LSPAD00, data = spurge),
  lm(dspad0t1 ~ LSPAD00 +LLAC01, data = spurge),
  LNIG00 %~~% LLAC00,
  spurge
  )

coefs(spurge_mod)

# Read Arkema data
arkema <- read.csv("../Data/arkema.csv")
arkema$site_trans <- paste(arkema$site, arkema$transect, sep="_")

#plot it
qplot(year, percent.algae, data=arkema, color = factor(transect)) +
  facet_wrap(~site) +
  geom_line()

#Model
algae_mod <- lm(percent.algae ~ frond.density, data = arkema)
invert_mod <- lm(percent.inverts ~ percent.algae, data = arkema)

arkema.sem <- psem(
  algae_mod,
  invert_mod,
  arkema
)


#Look at residuals by site
res_df <- data.frame(site_trans=arkema$site_trans,
                     res = residuals(algae_mod))

par(mfrow=c(2,2))
for(asite in levels(res_df$site_trans)){
  subdata <- subset(res_df, res_df$site_trans==asite)
  acf(subdata$res, cex.lab=1.3, main=asite)
}
par(mfrow=c(1,1))

summary(arkema.sem)

# Incorporate temporal autocorrelation

algae_mod_ac <- gls(percent.algae ~ frond.density, data = arkema,
                    correlation = corAR1(form = ~ year | site_trans))

invert_mod_ac <- gls(percent.inverts ~ percent.algae, data = arkema,
                     correlation = corAR1(form = ~ year | site_trans))

arkema.sem2 <- psem(
  algae_mod_ac,
  invert_mod_ac,
  arkema
)

coefs(arkema.sem2)

# Add site/transect grouping variable
arkema.sem3 <- psem(
  lme(percent.algae ~ frond.density, random = ~ 1 | site/transect,
      correlation = corCAR1(form = ~ year | site/transect),
      data =  arkema),
  
  lme(percent.inverts ~ percent.algae, random = ~ 1 | site/transect,
      correlation = corCAR1(form = ~ year | site/transect),
      data = arkema),
  
  data = arkema
)


coefs(arkema.sem3)
dSep(arkema.sem3)

summary(arkema.sem3)

AIC(arkema.sem2, arkema.sem3)

# Look at conflation of
ggplot(arkema, aes(x = site, y = mid.canopy)) +
  geom_bar(stat = "identity") +
  theme_bw(base_size = 24) +
  theme(legend.position = "none")

### Mixed Model approach