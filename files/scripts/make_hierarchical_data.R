library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(31415)

toplevel <- data.frame(site = letters[1:10], 
                       winter_temp = rnorm(10, 3),
                       winter_upwelling = rnorm(10)) %>%
  mutate(winter_phytoplankton = rnorm(10, winter_upwelling + winter_temp*1.4+4))


lowerlevel <- crossing(toplevel, s = 1:10) %>%
  select(-s) %>%
  mutate(algal_cover = abs(rnorm(100, 50, 25))) %>%
  group_by(site) %>%
  mutate(site_int = rnorm(1,0,5),
         ostracod_abund = rnorm(10, site_int + 10 + 0.3*algal_cover + 
                                  4*winter_temp + 5*winter_phytoplankton, 5)) %>%
  select(-site_int) %>%
  ungroup()

#take a look
qplot(winter_temp, winter_phytoplankton, color=site, data = lowerlevel)
qplot(algal_cover, ostracod_abund, color=site, data = lowerlevel)
qplot(winter_phytoplankton, ostracod_abund, color=site, data = lowerlevel)
qplot(winter_temp, ostracod_abund, color=site, data = lowerlevel)


library(lme4)
library(lmerTest)
a <- lmer(ostracod_abund ~ winter_temp + 
            winter_phytoplankton + algal_cover + (1|site), 
    data = lowerlevel)

summary(a)


write.csv(lowerlevel, "../data/ostracod_plotlevel.csv", row.names=FALSE)
write.csv(toplevel, "../data/ostracod_sitelevel.csv", row.names=FALSE)
