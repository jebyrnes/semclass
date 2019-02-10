#------------------
#' Bayesian SEM
#'
#' @author Jarrett Byrnes
#' 
#'@date 2018-03-15
#'
#------------------

library(brms)
library(piecewiseSEM)
library(ggplot2)

#Load the data
keeley <- read.csv("../data/Keeley_rawdata_select4.csv")

#The models
rich_mod <- bf(rich ~ firesev + cover)
cover_mod <- bf(cover ~ firesev)


#wrap it together in one model
k_fit_brms <- brm(rich_mod +
                    cover_mod + 
                    set_rescor(FALSE), 
                  data=keeley,
                  cores=4, chains = 2)
#check it
plot(k_fit_brms)

pp_check(k_fit_brms, resp="rich") +
  scale_color_manual(values=c("red", "black"))

# get the coefs
pp_check(k_fit_brms, resp="cover") +
  scale_color_manual(values=c("red", "black"))

# get the coefs
summary(k_fit_brms)

#plot to check convergence
plot(k_fit_brms)

#show that WAIC is the sum of pieces

rich_fit <- brm(rich_mod,
                data=keeley,
                cores=2, chains = 2)
cover_fit <- brm(cover_mod,
                 data=keeley,
                 cores=2, chains = 2)

WAIC(k_fit_brms)

WAIC(rich_fit)
WAIC(cover_fit)


# fit the fully mediated model
rich_mod_fullmed <- bf(rich ~  cover)

fit_brms_fullmed <- brm(rich_mod_fullmed +
                          cover_mod + 
                          set_rescor(FALSE), 
                        data=keeley,
                        cores=4, chains = 2)

#compare models
WAIC(k_fit_brms, fit_brms_fullmed)
WAIC(k_fit_brms, fit_brms_fullmed)

#show prediction and error propogation

#First, make new data
follow_this <- data.frame(firesev = 5)

#Get fitted sims of new data
# (predict gives you prediction intervals)
cover_fit <- fitted(k_fit_brms, newdata=follow_this,
                     resp = "cover", nsamples = 1000, 
                     summary = FALSE)

#make a new data frame for next variable, propogating simulation
#error into the predictor
follow_this_fit <- expand.grid(firesev = follow_this$firesev, 
                           cover = as.vector(cover_fit))

#second fit
rich_fit_values <- fitted(k_fit_brms, newdata=follow_this_fit,
                     resp = "rich", nsamples = 1000, 
                     summary = FALSE)

#remove excess simulations (from a 1000 x 1000 matrix)
rich_fit_values <- diag(rich_fit_values)

#posterior predictions
median(rich_fit_values)
posterior_interval(as.matrix(rich_fit_values))

ggplot(data = as.data.frame(rich_fit_values), aes(x=rich_fit_values)) + 
  geom_density(fill="blue", alpha=0.3) +
  theme_bw(base_size=17)

#now show prediction
cover_pred <- predict(k_fit_brms, newdata=follow_this,
                     resp = "cover", nsamples = 1000, 
                     summary = FALSE)

follow_this_pred <- expand.grid(firesev = follow_this$firesev, 
                           cover = as.vector(cover_pred))

rich_pred <- predict(k_fit_brms, newdata=follow_this_pred,
                    resp = "rich", nsamples = 1000, 
                    summary = FALSE)

rich_pred <- as.matrix(diag(rich_pred))


median(rich_pred)
posterior_interval(rich_pred)

########
#prediction curves
########

newdata_curve = data.frame(firesev = c(1,10))

#fitted
cover_pred_curve <- fitted(k_fit_brms, newdata=newdata_curve,
                     resp = "cover", nsamples = 1000, 
                     summary = FALSE)

#with error
cover_pred_error_curve <- predict(k_fit_brms, newdata=newdata_curve,
                     resp = "cover", nsamples = 1000, 
                     summary = FALSE)

plotdata <- data.frame(firesev_start=1,
                       firesev_end = 10,
                       fit_start = cover_pred_curve[,1],
                       fit_end = cover_pred_curve[,2],
                       pred_start = cover_pred_error_curve[,1],
                       pred_end = cover_pred_error_curve[,2]
                       )

plotdata_median <- data.frame(firesev_start=1,
                              firesev_end = 10,
                              pred_start = median(cover_pred_curve[,1]),
                              pred_end = median(cover_pred_curve[,2]))

ggplot(plotdata, mapping=aes(x=firesev_start, xend = firesev_end,
                             y = pred_start, yend = pred_end)) +
  geom_segment(color="black",
               alpha = 0.1) +
  geom_segment(mapping=aes(y = fit_start, yend = fit_end), 
               color="lightblue",
               alpha = 0.1) +
  geom_segment(data = plotdata_median, color="red", lwd=1.5) +
  theme_bw(base_size=17) +
  xlab("firesev") + ylab("cover")


#prediction of second variable
newdata2 <- expand.grid(firesev = newdata_curve$firesev, cover = as.vector(cover_pred))

newdata2_error <- expand.grid(firesev = newdata_curve$firesev, cover = as.vector(cover_pred_error_curve))

rich_pred <- fitted(k_fit_brms, newdata=newdata2,
                    resp = "rich", nsamples = 1000, 
                    summary = FALSE)

rich_pred_error <- fitted(k_fit_brms, newdata=newdata2_error,
                    resp = "rich", nsamples = 1000, 
                    summary = FALSE)


#to minimize excess uncertainty
rich_pred <- as.matrix(diag(rich_pred))
rich_pred_error <- as.matrix(diag(rich_pred_error))

#visualize
par(mfrow=c(1,2))
plot(density(as.vector(rich_pred)), main = "fitted",
     xlim = c(25, 65))
plot(density(as.vector(rich_pred_error)), main = "prediction",
     xlim = c(25, 65))
par(mfrow=c(1,1))
