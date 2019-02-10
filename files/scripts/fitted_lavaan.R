library(lavaan)


predict_lavaan <- function(fit, newdata = NULL, 
                           simCI=FALSE, nsims=100, ci=0.95,
                           interval="confidence"){
  
  stopifnot(inherits(fit, "lavaan"))
  
  #Make sure we can use this
  if(!inspect(fit, "meanstructure")) stop("Need to supply meanstructure = TRUE in fit\n")
  if(is.null(newdata)){
    newdata <- data.frame(inspect(fit, "data"))
    names(newdata) <- lavNames(fit)
  }
  
  if(length(lavNames(fit, type="lv"))!=0) stop("Does not currently work with latent variables\n")
  
  #check for new data
  if(sum(!(lavNames(fit, type="ov.x") %in% names(newdata)))>0) stop("Not all exogenous variables supplied!")
  
  #Add some new columns to newdata
  newdata$Intercept <- 1
  newdata[lavNames(fit, "ov.nox")] <- 0
  
  #get coefficients
  mod_df <- make_mod_df(fit)
  
  #Drop covariances or derived quantities
  mod_df <- mod_df[-which(mod_df$op=="~~"),]
  if(length(which(mod_df$op==":="))>0) mod_df <- mod_df[-which(mod_df$op==":="),]
  mod_df[which(mod_df$op=="~1"),]$rhs <- "Intercept"
  
  #get rid of exogenous on lhs
  mod_df <- mod_df[-which(mod_df$exo==1),]
  
  #Order by lhs
  mod_df <- mod_df[sort(mod_df$lhs, index.return=TRUE)$ix,]
  
  #let us know which variables on the rhs are exogenous
  mod_df$ord <- 0
  mod_df[which(!(mod_df$rhs %in% mod_df$lhs)),]$ord <- 1
  
  #Make a "order"
  ord_current <- 1
  while(sum(mod_df$ord==0)>0){
    for(r in unique(mod_df$lhs)){
      val <-  sum(mod_df[which(mod_df$lhs==r),]$ord==0)
      if(val==0) {
        mod_df[which(mod_df$lhs==r),]$ord <- ord_current
        
        if(sum(mod_df$rhs==r)>0)
          mod_df[which(mod_df$rhs==r),]$ord <- ord_current+1
      }
    }
    ord_current <- ord_current +1
  }
  
  #correct for ragged ordering
  for(r in unique(mod_df$lhs)){
    mod_df[which(mod_df$lhs==r),]$ord <- max(mod_df[which(mod_df$lhs==r),]$ord)
  }
  
  #sort by order 
  mod_df <- mod_df[sort(mod_df$ord, index.return=TRUE)$ix,]
  
  
  #get the mean estimates
  ret <- get_predicted_values_lavaan(mod_df, newdata)
  if(simCI) ret <- cbind(ret, get_sim_ci_lavaan(fit, mod_df, newdata, nsims, ci, interval))
  
  
  return(as.data.frame(apply(ret, 2, as.numeric)))
  
}

make_mod_df <- function(fit){
  
 data.frame(lhs = fit@ParTable$lhs,
            op = fit@ParTable$op,
            rhs = fit@ParTable$rhs,
            exo = fit@ParTable$exo,
            est = fit@ParTable$est,
            se = fit@ParTable$se,
            label=fit@ParTable$label,
            stringsAsFactors=FALSE)
}

get_predicted_values_lavaan <- function(mod_df, newdata, err=FALSE){
  #now do the fitting in order
  fit_df <- data.frame(base = rep(1, nrow(newdata)))
  
  for(r in unique(mod_df$lhs)){
    subdf <- subset(mod_df, mod_df$lhs==r)
    if(err){
      err_est <-  subdf[which(subdf$op == "~~"),]$est
      subdf <- subdf[-which(subdf$op == "~~"),]
    }
    
    #make a formula
    rhs <- paste0(subdf$rhs, collapse=" + ")
    form <- as.formula(paste0(r, " ~ ", rhs))
    
    #use formula to get right part of the data in right format
    mod_mat <- model.matrix(form, newdata)[,-1]
    new_val <- mod_mat %*% subdf$est
    
    #in case we are doing prediction error
    if(err) new_val <- rnorm(length(new_val), new_val, sd=sqrt(err_est))
    
    fit_df[[r]] <- new_val
    newdata[[r]] <- new_val
  }
  
  return(as.vector(fit_df[,-1]))
}

get_sim_ci_lavaan <- function(fit, mod_df, newdata, nsims, ci, interval="confidence"){
  #in case we're doing something with SEs
  eqns <- paste0(mod_df$lhs, mod_df$op, mod_df$rhs)
  eqns <- gsub("Intercept", "", eqns)
  
  #If there are labels, we'll need to include them in 
  #eqns for proper vcov matrix ordering
  if(sum(mod_df$label != "")>0){
    eqns[which(mod_df$label != "")] <- mod_df$label[which(mod_df$label != "")]
  }
  
  
  #get the proper ordered variance-covariance matrix
  vmat <- vcov(fit)
  vmat <- vmat[-grep("~~", rownames(vmat)),-grep("~~", colnames(vmat))]
  
  
  vcov_mat <- vmat[which(rownames(vmat) %in% eqns),]
  vcov_mat <- vcov_mat[,which(colnames(vcov_mat) %in% eqns)]
  
  vcov_mat <- vcov_mat[match(eqns, rownames(vcov_mat)),]
  vcov_mat <- vcov_mat[,match(eqns, colnames(vcov_mat))]
  
  #build in the residuals
  if(interval=="prediction"){
    err_vmat <-  vcov(fit)
    err_vmat <- err_vmat[grep("~~", rownames(err_vmat)),]
    err_vmat <- err_vmat[,match(colnames(vcov_mat), colnames(err_vmat))]

    err_vcmat <-  vcov(fit)
    err_vcmat <- err_vcmat[grep("~~", rownames(err_vcmat)),
                          grep("~~", colnames(err_vcmat))]
    
    vcov_mat <- rbind(vcov_mat,err_vmat)
    vcov_mat <- cbind(vcov_mat, t(cbind(err_vmat, err_vcmat)))
    
    #add estimates to mod_df
    err_terms <- make_mod_df(fit)
    err_terms <- err_terms[which(err_terms$rhs == err_terms$lhs),]
    err_terms <- err_terms[which(err_terms$exo == 0),]
    err_terms$ord <- 0
    mod_df <- rbind(mod_df, err_terms)
    
    #get rid of correlated error terms
    corr_err <- with(mod_df, which(op=="~~" & lhs != rhs))
    if(length(corr_err)>0) mod_df <- mod_df[-corr_err,]
    
    #reorder vcov mat
    eqns <- paste0(mod_df$lhs, mod_df$op, mod_df$rhs)
    eqns <- gsub("Intercept", "", eqns)
    vcov_mat <- vcov_mat[match(eqns, rownames(vcov_mat)),]
    vcov_mat <- vcov_mat[,match(eqns, colnames(vcov_mat))]
  }
  
  #If we are getting simulated CIs
  est_sim <- mvtnorm::rmvnorm(nsims, mod_df$est, sigma=vcov_mat)
  sims <- lapply(1:nrow(est_sim), function(i){
    arow <- est_sim[i,]
    mdf <- mod_df
    mdf$est = arow
    sim = get_predicted_values_lavaan(mdf, newdata, err=(interval=="prediction"))
    sim$pt <- 1:nrow(sim)
    sim
  })
  
  sims <- do.call(rbind, sims)
  
  #seriously this is easier with dplyr,
  #but I'm executing package discipline
  
  lower_ci <- stats::aggregate(. ~ pt, data=sims, quantile, probs=(1-ci)/2)[,-1]
  names(lower_ci) <- paste0(names(sims)[-ncol(sims)], "_lower")
  
  upper_ci <- stats::aggregate(. ~ pt, data=sims, quantile, probs=1-(1-ci)/2)[,-1]
  names(upper_ci) <- paste0(names(sims)[-ncol(sims)], "_upper")
  
  
  return(cbind(lower_ci, upper_ci))
}


fitted_lavaan <- function(fit){
  predict_lavaan(fit)
}

residuals_lavaan <- function(fit){
  fitted_vals <- fitted_lavaan(fit)
  
  rawdata <- data.frame(inspect(fit, "data"))
  names(rawdata) <- lavNames(fit)
  
  res <- data.frame(base = rep(1, nrow(rawdata)))
  for(vals in names(fitted_vals)){
    res[[vals]] <- rawdata[[vals]] - fitted_vals[[vals]] 
  }
  
  return(res[,-1])
}


############################################
########### EXAMPLE CODE
############################################
# 
# iris_mod <- "
# Sepal.Width ~ Sepal.Length
# Petal.Width ~ Sepal.Width + Petal.Length
# 
# "

# iris_fit <- sem(iris_mod, data=iris, meanstructure = TRUE)

# fitted_lavaan(iris_fit)
# 
# res <- residuals_lavaan(iris_fit)
# 
# par(mfrow=c(1,2))
# for(i in 1:2) {
#   qqnorm(res[,i], main=names(res)[i])
#   qqline(res[,i])
# }
# par(mfrow=c(1,1))
# # 
# newdata <- data.frame(expand.grid(Sepal.Length=1:10,
#                        Petal.Length=c(1,10)))
# pred <- predict_lavaan(iris_fit, newdata=newdata)
# 
# plot(Petal.Width ~ Sepal.Length, data=cbind(newdata, pred))
# 
# # library(mvnormtest)
# # mshapiro.test(t(res))
# #
# pred2 <- predict_lavaan(iris_fit, newdata=newdata, simCI=TRUE)
# pred3 <- predict_lavaan(iris_fit, newdata=newdata, simCI=TRUE, interval="prediction")
# pred2 <- cbind(newdata, pred2)
# pred3 <- cbind(newdata, pred3)
# 
# library(ggplot2)
# ggplot(pred2, mapping=aes(x=Sepal.Length, y=Petal.Width,
#                           ymin = Petal.Width_lower,
#                           ymax=Petal.Width_upper)) +
#   geom_ribbon(alpha=0.4, fill="blue") +
#   geom_ribbon(alpha=0.4, fill="grey", data=pred3) +
#   geom_line() +
#   facet_wrap(~Sepal.Length, scale="free_y") +
#   theme_bw()
# 
# 
# #show fit v. prediction CI
# ggplot(pred2, mapping=aes(x=Sepal.Length, y=Petal.Length,
#                           size = Petal.Width_upper - Petal.Width_lower,
#                           color = Petal.Width)) +
#   geom_point(alpha=0.4) +
#   geom_point(alpha=0.4, color="black", data=pred3, shape=1) +
#   theme_bw() +
#   scale_color_gradient(low="blue", high="red")


# plot(pred2$Petal.Width_lower, pred3$Petal.Width_lower)
