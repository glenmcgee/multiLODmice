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
#'
#' @seealso \code{\link{pool_ests}} and \code{\link{pool_lm}}
#'
#' @examples
#' set.seed(123)
#' require(crch)
#' require(censReg)
#' require(mvtnorm)
#'
#' ## generate complete data
#' n <- 5000
#' dat <- gendat_multiLOD(n,nonlinear=TRUE)
#' X <- dat$X
#' df <- dat$df
#' LODs <- dat$LODs
#'
#' ## model fit with no censoring
#' g_gold <- lm(df$y~X1+X2+X3+X4,data=X) ## gold standard
#'
#' ## complete cases only
#' g_cc <- lm(y~X1+X2+X3+X4,data=df)
#'
#' ## naive imputation: LOD/sqrt(2)
#' df_naive <- df
#' df_naive$X1[!is.na(LODs$X1)] <- LODs$X1[!is.na(LODs$X1)]/sqrt(2)
#' df_naive$X2[!is.na(LODs$X2)] <- LODs$X2[!is.na(LODs$X2)]/sqrt(2)
#' g_naive <- lm(y~X1+X2+X3+X4,data=df_naive)
#'
#' ## MI
#' kk <- 20
#' df_imp <- multiLODmice(data = df, data.lod = LODs,mi.m = kk)
#' g_imp <- pool_lm(y~X1+X2+X3+X4,df_imp)
#' g_imp <- pool_ests(g_imp)
#'
#' ## compare percent bias
#' beta <- c(0,1,1,1,1)
#' round(100*(coef(g_gold)-beta))
#' round(100*(coef(g_cc)-beta))
#' round(100*(coef(g_naive)-beta))
#' round(100*(g_imp$ests-beta))
#'
#' ## compare uncertainty
#' round(sqrt(diag(vcov(g_gold))),3)
#' round(sqrt(diag(vcov(g_cc))),3)
#' round(sqrt(diag(vcov(g_naive))),3)
#' round(sqrt(diag(g_imp$vcov)),3)
#'
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
#' @return List of estimates vector and variance-covariance matrix.
#'
#' @seealso \code{\link{pool_ests}} and \code{\link{pool_gam}}
#'
#' @examples
#' set.seed(123)
#' require(crch)
#' require(censReg)
#' require(mvtnorm)
#'
#' ## generate complete data
#' n <- 5000
#' dat <- gendat_multiLOD(n)
#' X <- dat$X
#' df <- dat$df
#' LODs <- dat$LODs
#'
#' ## model fit with no censoring
#' g_gold <- lm(df$y~X1+X2+X3+X4,data=X) ## gold standard
#'
#' ## complete cases only
#' g_cc <- lm(y~X1+X2+X3+X4,data=df)
#'
#' ## naive imputation: LOD/sqrt(2)
#' df_naive <- df
#' df_naive$X1[!is.na(LODs$X1)] <- LODs$X1[!is.na(LODs$X1)]/sqrt(2)
#' df_naive$X2[!is.na(LODs$X2)] <- LODs$X2[!is.na(LODs$X2)]/sqrt(2)
#' g_naive <- lm(y~X1+X2+X3+X4,data=df_naive)
#'
#' ## MI
#' kk <- 20
#' df_imp <- multiLODmice(data = df, data.lod = LODs,mi.m = kk)
#' g_imp <- pool_lm(y~X1+X2+X3+X4,df_imp)
#' g_imp <- pool_ests(g_imp)
#'
#' ## compare percent bias
#' beta <- c(0,1,1,1,1)
#' round(100*(coef(g_gold)-beta))
#' round(100*(coef(g_cc)-beta))
#' round(100*(coef(g_naive)-beta))
#' round(100*(g_imp$ests-beta))
#'
#' ## compare uncertainty
#' round(sqrt(diag(vcov(g_gold))),3)
#' round(sqrt(diag(vcov(g_cc))),3)
#' round(sqrt(diag(vcov(g_naive))),3)
#' round(sqrt(diag(g_imp$vcov)),3)
#'
#'
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
#'
#' @seealso \code{\link{pool_lm}}
#'
#' @example
#' ## See \code{\link{pool_lm}}
#'
#' @export
pool_ests <- function(obj){
  mm <- length(obj) ## number of imputations

  ## get estimates and vcov estimates from fitted models
  ests <- t(sapply(obj,function(x) x$ests)) ## get matrix of all estimates
  vcovs <- lapply(obj,function(x) x$vcovs) ## get list of all vcov matrices

  ## rubin's rules
  MI_ests <- apply(ests,2,mean)
  MI_vcov <- Reduce('+', vcovs)/mm + (1+1/mm)*cov(ests) # average vcov + sample cov of estimates

  return(list(ests=MI_ests,
              vcov=MI_vcov))
}

