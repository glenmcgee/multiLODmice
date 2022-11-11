## pool_funs.R
## functions to pool estimates


#' Fit GAM to each imputed dataset and combine.
#'
#' Uses Rubin's rules to combine GAM fits across imputed datasets.
#' @param form model formula as in mgcv
#' @param impdat list of imputed datasets
#' @param fixedknots list of fixed knots to be used by all model fits
#' @param fam family (default="gaussian") to be passed to mgcv::gam()
#' @param combinemethod method for combining results. Default="betas" combines coefficients as in lm(). Otherwise combines fitted values
#'
#'
#' @return List of estimates vector and variance-covariance matrix.
#' @export
pool_gam <- function(form,impdat,fixedknots,fam="gaussian",combinemethod="betas"){
  ## get formula for lm
  form <- as.formula(form)

  ## fit gam to each imputed dataset
  fits <- lapply(impdat,function(dat) mgcv::gam(form,data=dat,
                                                knots=fixedknots,
                                                family=fam,
                                                method="REML")) # return list of gam fits

  ## return list of ests and vcovs for each model fit
  if(combinemethod=="betas"){ ## apply rubins rules to underlying coefficients
    res <- lapply(fits,function(fit) list(ests=coef(fit),
                                          vcovs=vcov(fit)))
  }else{ ## apply rubins rules to fitted vals/curves directly

  }

}

#' Fit LM/GLM to each imputed dataset and combine.
#'
#' Uses Rubin's rules to combine LM/GLM fits across imputed datasets.
#' @param form model formula as in lm/glm
#' @param impdat list of imputed datasets
#' @param fam family (default="gaussian") to be passed to glm()
#'
#'
#' @return List of estimates vector and variance-covariance matrix.
#' @export
pool_lm <- function(form,impdat,fam="gaussian"){

  ## get formula for lm
  form <- as.formula(form)

  ## fit model to each imputed dataset
  if(fam=="gaussian"){
    fits <- lapply(impdat,function(dat) lm(form,data=dat)) # return list of model fits
  }else{
    fits <- lapply(impdat,function(dat) glm(form,data=dat,family=fam)) # return list of model fits
  }

  ## return list of ests and vcovs for each model fit
  res <- lapply(fits,function(fit) list(ests=coef(fit),
                                        vcovs=vcov(fit)))
}

#' Pool generic estimates.
#'
#' Internal function for pooling generic estimates via Rubin's rules.
#' @param obj model formula as in mgcv
#'
#'
#' @return List of estimates vector and variance-covariance matrix.
#' @export
pool_ests <- function(obj){
  mm <- length(obj) ## number of imputations

  ## get estimates and vcov estimates from fitted models
  ests <- t(sapply(obj,function(x) x$ests)) ## get matrix of all estimates
  vcovs <- lapply(obj,function(x) x$vcovs) ## get list of all vcov matrices

  ## rubin's rules
  MI_ests <- apply(ests,2,mean)
  MI_vcov <- Reduce('+', vcovs)/mm + cov(ests) # average vcov + sample cov of estimates

  return(list(ests=MI_ests,
              vcov=MI_vcov))
}

