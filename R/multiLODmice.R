### MODIFIED MICE FUNCTIONS FOR LEFT-CENSORING WITH MULTIPLE LODs
## comes from doMIsaul package by Faucheux et al (https://github.com/LilithF/doMIsaul)
## method is based on Lapidus et al (2014; Statistics in Medicine), extended to multiple/varying LODs
## GNU General Public License v3.0 (https://github.com/LilithF/doMIsaul/blob/main/LICENSE.md)



#' Base function for imputing multiply left censored data with MICE,
#'
#' Based on \code{Lapidus et al.} extended to allow multiple LODs. Modified version of .cens.draw3 function from doMIsaul package.
#'
#' @keywords internal
#'
#' @param y Vector to be imputed
#' @param ry Logical vector of length \code{length(y)} indicating the the subset
#'   \code{y[ry]} of elements in y to which the imputation model is fitted.
#'   The \code{ry} generally distinguishes the observed (\code{TRUE}) and
#'   missing values (\code{FALSE}) in y.
#' @param x Numeric design matrix with \code{length(y)} rows with predictors
#'   for \code{y}. Matrix \code{x} may have no missing values.
#' @param lod.j Vector of varying LODs (left censoring values). -Inf for non censored values.
#' @param ... Other named arguments.
#'
#' @return parameters
.multicens.draw3 <- function(y,  ## passed yy from above: all observed OR censored values
                             ry, ## passed ryy from above:  indicator of it being observed OR censored
                             x,  ## these are all the other data (excluding the LOD columns)
                             lod.j, ## vector of lod.js (and -Inf for observed vals & non-censored missing vals)# single LOD
                             ...) {

  x <- as.data.frame(x)
  ## subset to observed vals AND censored values; exclude truly missing values
  x <- x[ry, ]
  y <- y[ry]
  lod.j <- lod.j[ry]

  ## formula for imputation model (regress var on all others)
  form <- as.formula(paste("outcome_y", paste(names(x)[-1], collapse = " + "), sep = " ~ "))

  df <- cbind(outcome_y=y,x)

  ## fit imputation model using crch
  fit3 <- crch::crch(form,data=df, dist = "gaussian",left=lod.j,truncated=FALSE) ## truncated=FALSE gives correct CENSORED regression
  # fit3 <- censReg::censReg(form, left = lod.j, data = x) ## censReg only allows single lod.j value


  ## draw beta ans logsigma from posterior/sampling distn
  draw <- MASS::mvrnorm(1, mu = coef(fit3), Sigma = vcov(fit3))
  # draw <- MASS::mvrnorm(1, mu = fit3$estimate, Sigma = solve(-fit3$hessian)) # censReg version

  ## save draw of beta as well as sigma for imputations
  parm <- list(beta = as.numeric(draw[-length(draw)]),
               sigma = exp(as.numeric(draw[length(draw)])))

  return(parm)
}

#' Impute left censored data with MICE
#'
#' Function based on \code{Lapidus et al.} for imputing left-censored data with mice; extended to varying LODs. Modified version of mice.impute.cens function from doMIsaul package.
#' @param y Vector to be imputed
#' @param ry Logical vector of length \code{length(y)} indicating the the subset
#'   \code{y[ry]} of elements in y to which the imputation model is fitted.
#'   The \code{ry} generally distinguishes the observed (\code{TRUE}) and
#'   missing values (\code{FALSE}) in y.
#' @param x Numeric design matrix with \code{length(y)} rows with predictors
#'   for \code{y}. Matrix \code{x} may have no missing values.
#' @param lod.j Vector of varying LODs (left censoring values). -Inf for non censored values.
#' @param wy Logical vector of length \code{length(y)}. A \code{TRUE} value
#'   indicates locations in \code{y} for which imputations are created.
#' @param ... Other named arguments.
#'
#' @return Vector with imputed data, same type as \code{y}, and of length
#'   \code{sum(wy)}.
#' @export
mice.impute.multicens <- function(y, ## the outcome for this imputation model
                             ry, ## indicator of being observed or NA (censored values are considered NA so that they can be imputed)
                             x, ## all other covariates, LOD columns included
                             lod.j, ## vector of LODs
                             excludenames, ## vector of names to be excluded
                             wy = NULL, ...) {

  ## define wy to be an indicator of observations to be imputed
  if (is.null(wy)){
    wy <- !ry ## set to the observations that are NA (i.e. not observed)---this includes censored values
  }

  ## x is just the predictors
  x <- cbind(1,as.matrix(x))
  x <- x[,!(colnames(x) %in% excludenames)]

  ## define yy to be the observed variable AND ALSO the LOD where it has been censored
  yy <- y
  yy[is.finite(lod.j)] <- lod.j[is.finite(lod.j)]

  ## define ryy to be an indicator of observing the variable OR observing the LOD (i.e. FALSE ONLY IF TRULY MISSING)
  ryy <- ry ## indicator that values are fully observed
  ryy[is.finite(lod.j)] <- TRUE # EDIT: now just set to TRUE when LOD is finite ## treat censored values as observed for fitting imputation model

  ## define method to fit left censored imputation model and get posterior draw of betas+sigma^2 from imputation model
  doMI.multicens.draw3 <- get(".multicens.draw3", envir = asNamespace("multiLODmice"), # use envir = as.environment(1) to instead search the global environment
                         inherits = FALSE)

  ## call function to fit l-cens imputation model
  ### this function uses yy and ryy, meaning it INCLUDES censored values, and respects the censoring
  parm <- doMI.multicens.draw3(y = yy, ry = ryy, x = x, lod.j = lod.j) # EDIT: pass entire vector xx as LODs

  ## AMONG the observations to be imputed, which ones are really censored (and hence equal to lod.j)
  wyy <- wy[wy]
  wyy[!is.finite(lod.j[wy])] <- FALSE # EDIT: ## set to FALSE if not actually a censored val (hence FALSE values should be truly missing)

  ## for all values to be imputed, impute from Normal distribution using beta and sigma from above
  .draw <- as.numeric (x[wy, ] %*% parm$beta +
                         rnorm(sum(wy)) * parm$sigma)

  ## for values to be imputed that were really censored, draw from a TRUNCATED normal using beta and sigma from above
  .draw[wyy] <- truncnorm::rtruncnorm(
    sum(wyy), a = -Inf, b = lod.j[wy & is.finite(lod.j) ], #EDIT: lod.j now a vector
    mean = x[wy & is.finite(lod.j), ] %*% parm$beta, #EDIT: # mean = x[wy & xx == lod.j, ] %*% parm$beta,
    sd = parm$sigma
  )

  .draw <- as.matrix(.draw)
  return(.draw)
}
environment(mice.impute.multicens) <- environment(mice::mice.impute.norm)



#' Main function for mice with multiply left censored data.
#'
#' Calls mice with custom method for multiply left censored data. Based on \code{Lapidus et al.} for imputing left-censored data with mice; extended to varying LODs. Modified version of MImpute_lcens function from doMIsaul package.
#'
#' @param data dataset including some missing or censored values
#' @param data.lod dataset containing columns corresponding to LOD values (left censors) and NA elsewhere (for non-censored observations or truly missing values). Column names must match columns from \code{data} with censored values.
#' @param mi.m Number of imputed datasets
#' @param excludeLODimp vector of names of variables to be excluded from truncated LOD imputation model
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
#' ## impute data
#' kk <- 20
#' df_imp <- multiLODmice(data = df, data.lod = LODs,mi.m = kk)
#'
#' ## check that imputed values are correctly below the LOD
#' sum(sapply(df_imp,function(dat) sum(dat$X1[is.na(df$X1) & seq(1,n)<(n/2)]> (-2))))==0 ## first LOD for X1
#' sum(sapply(df_imp,function(dat) sum(dat$X1[is.na(df$X1) & seq(1,n)>(n/2)]> (-1))))==0 ## second LOD for X1
#' sum(sapply(df_imp,function(dat) sum(dat$X2[is.na(df$X2) & seq(1,n)<(n/2)]> (-1))))==0 ## first LOD for X2
#' sum(sapply(df_imp,function(dat) sum(dat$X2[is.na(df$X2) & seq(1,n)>(n/2)]> (-2))))==0 ## second LOD for X2
#'
#' ## compare imputations to true values
#' mean_imp1 <- apply(sapply(df_imp,function(dat) dat$X1[is.na(df$X1)]),1,mean)
#' mean_imp2 <- apply(sapply(df_imp,function(dat) dat$X2[is.na(df$X2)]),1,mean)
#'
#' ## check accuracy
#' summary(lm(X[is.na(df$X1),1]~mean_imp1))
#' summary(lm(X[is.na(df$X2),2]~mean_imp2))
#'
#'
#' @export
multiLODmice <- function(data,      ## main dataset
                         data.lod,  ## EDITED: LOD columns indicating LOD (for censored values) or -Inf (for un-censored values--these are either observed or truly missing) ## should have same names as corresponding data columns
                         mi.m,      ## number of imputed datasets
                         excludeLODimp=NULL, ## vector of names of variables to be excluded from truncated LOD imputation model
                         maxit = 10,
                         return.midsObject = FALSE){

  data.lod[is.na(data.lod)] <- -Inf ## convert NAs in data.lod to an LOD of -Inf

  ## select mice method for imputation
  my.method <- sapply(colnames(data), function(i){
    if (sum(is.na(data[, i])) > 0){
      if (is.element(i, colnames(data.lod))){
        "multicens" ## if there are NAs AND there is an LOD column, assign lcens mice method for left censoring
      } else {
        "pmm"## if there are NAs and NO LOD column, assign default imputation method
      }
    } else {
      "" ## if there are no NAs dont do anything
    }

  }, USE.NAMES = TRUE)

  ## pass entire vector of time-varying LODs from data.lod
  my.blots <- sapply(names(my.method), function(i){
    if (my.method[i] == "multicens") {
      list(lod.j = data.lod[, i], ## passing the entire vector of LODs (including -Inf for non-censored values)
           excludenames=excludeLODimp)
    } else {
      list()
    }
  }, USE.NAMES = TRUE, simplify = FALSE)

  ## call mice using custom methods above
  res <- mice::mice(data, m = mi.m, method = my.method,
                    maxit = maxit, blots = my.blots, print = FALSE)

  imp.ret <- mice::complete(res, "all")

  if(return.midsObject){
    imp2 <- mice::complete(res, "long", include = TRUE)
    if (is.numeric(mice.log)){
      imp2[, -c(1:2)] <- exp(imp2[, -c(1:2)])
    }
    imp2 <- imp2[, 1:ncol(data)]
    MIDS <- mice::as.mids(imp2)

    return(list(
      imputed.data = imp.ret,
      mids.obj = MIDS))
  } else {
    return(imp.ret)
  }

}

#
# #### TESTING IT OUT:
# set.seed(1)
# toyraw <- iris[, 1:4]
# toy <- toyraw
# toy <- toy[sample(nrow(toy)),]
#
# Censored <- toy[,1:2]
# c1 <- rep(quantile(toy[,1],c(0.1,0.25,0.45)),each=50)
# Censored[toy[,1]<c1,1] <- c1[toy[,1]<c1]
# Censored[toy[,1]>=c1,1] <- NA
# c2 <- rep(quantile(toy[,2],c(0.1,0.25,0.45)),each=50)
# Censored[toy[,2]<c2,2] <- c2[toy[,2]<c2]
# Censored[toy[,2]>=c2,2] <- NA
#
# toy[toy[,1]<c1,1] <- NA
# toy[toy[,2]<c2,2] <- NA
#
#
# toy.imp <- multiLODmice(data = toy, data.lod = Censored,
#                          mi.m = 5)
#
#
# i=1
# j=3
# range=c((j-1)*50 + (1:50))
# max(c(toy.imp$`1`[range[which(is.na(toy[range,i]))],i],
#       toy.imp$`2`[range[which(is.na(toy[range,i]))],i],
#       toy.imp$`3`[range[which(is.na(toy[range,i]))],i],
#       toy.imp$`4`[range[which(is.na(toy[range,i]))],i],
#       toy.imp$`5`[range[which(is.na(toy[range,i]))],i]))
# quantile(toyraw[,i],c(0.1,0.25,0.45)[j])
# hist(c(toy.imp$`5`[range[which(is.na(toy[range,i]))],i]),xlim=range(toyraw[,i]))
# hist(c(toy.imp$`2`[range[which(!is.na(toy[range,i]))],i]),xlim=range(toyraw[,i]))
