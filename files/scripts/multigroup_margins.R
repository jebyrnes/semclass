# Standardized Coefficients
library(margins)

#For all
getMargins.psem <- function(x, at){
  x <- removeData(x, formula = 1)
  
  ret <- lapply(x,getMargins, at = at)
  
  #make a df
  ret <- do.call(rbind, ret)
  
  #remove intercepts
  ret[-which(ret$factor %in% names(at)),]
}

getMargins <- function(x, at = NULL) {
  ret <- summary(margins(x, at = at))
  
  ret <- cbind(response = rep(as.character(attr(terms(x), "variables")[2]), nrow(ret)), ret)
  ret[,1:4]
}


getSDVars <- function(object, by){
  require(tidyr)
  dat <- GetData(object)
  modList <- removeData(object, formula = 1)
  formulaList <- get.formula.list(modList)
  vars <- unique(unlist(sapply(formulaList, all.vars)))
  
  dat <- dat[,which(names(dat) %in% vars)]
  
  #get groups
  by_list <- lapply(by, function(x) dat[[x]])
  names(by_list) <- by
  
  #only numerics
  dat <- dat[,which(sapply(dat, class) %in% c("integer", "numeric"))]

  #deal with integer interactions
  bycols <- which(names(dat) %in% by)
  if(length(bycols>0)) dat <- dat[,-bycols]
  
  sds <- aggregate(dat, by = by_list, FUN = sd)
  
  sds %>% gather(variable, value, -!!(by))

}

getStdByGroup <- function(object, by){
  require(dplyr)
  dat <- GetData(object)
  
  at <- list(unique(dat[[by]]))
  names(at) <- by
  
  pcoefs <- getMargins.psem(object, at = at)
  sd_df <- getSDVars(object, by)
  
  left_join(pcoefs, sd_df, by = c("factor" = "variable", by)) %>%
    rename(sd_x = value) %>%
    left_join(sd_df, by = c("response" = "variable",  by)) %>%
    rename(sd_y = value) %>%
    
    #add std coefs
    mutate(std_coef = AME * sd_x/sd_y)
  
}