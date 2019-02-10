library(mvtnorm)
library(tidyverse)
library(ggplot2)

set.seed(2002)
a <- rmvnorm(100, c(0,0), sigma=matrix(c(1, 0.9, 0.9, 1), byrow=TRUE, ncol=2))
set.seed(2002)
b <- rmvnorm(100, c(0,0), sigma=matrix(c(1, 0.4, 0.4, 1), byrow=TRUE, ncol=2))


#get the slopes right
z <- lm(X2 ~ X1, data=data.frame(a))
z2 <- lm(X2 ~ X1, data=data.frame(b))
z
z2

b[,2] <- b[,2]  * coef(z)[2]/coef(z2)[2]
lm(X2 ~ X1, data=data.frame(a))
lm(X2 ~ X1, data=data.frame(b))


dat <- data.frame(rbind(a,b))
names(dat) <- c("x", "y")
dat <- mutate(dat, r=c(rep("r = 0.9", 100), rep("r = 0.4",100)))

qplot(x,y,data = dat) +
  stat_smooth(method="lm") +
  theme_bw(base_size=17) +
  facet_wrap(~r) +
  ggtitle(str_c("Slope = ", round(coef(z)[2],3)))


dat2 <- data.frame(rbind(b, b))    
names(dat2) <- c("x", "y")
dat2$grp <- c(rep(1,100), rep(2,100))
dat2[101:200,2] <- dat2[101:200,2]*3
dat2 <- dat2 %>%
  group_by(grp) %>%
  mutate(slope = coef(lm(y~x))[2]) %>%
  mutate(r = cor(x,y)) %>%
  mutate(slope = str_c("slope = ", round(slope,3)))

qplot(x,y,data = dat2) +
  stat_smooth(method="lm") +
  theme_bw(base_size=17) +
  facet_wrap(~slope) +
  ggtitle("r = 0.424")

