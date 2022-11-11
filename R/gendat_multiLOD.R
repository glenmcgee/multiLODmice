#' Generate toy data for testing functions.
#'
#' Function generates toy data with multiple LODs for testing functions.
#' @param n sample size
#' @param betas vector containing 5 coefficients
#' @param rho correlation for X's
#' @param LOD1 vector containing first and second LODs for X1
#' @param LOD2 vector containing first and second LODs for X2
#' @param nonlinear should the function use non-linear associations (default=FALSE)
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
#' head(dat$X)
#' head(dat$df)
#' head(dat$LODs)
#'
#' @return true data, observed data and LOD data
#' @export
gendat_multiLOD <- function(n=5000,
                            LOD1=c(-2,-1),LOD2=c(-1,-2),
                            rho=0.5,
                            betas=c(0,1,1,1,1),
                            nonlinear=FALSE){

  ## generate complete data
  X <- rmvnorm(n,
               mean=rep(0,4),
               sigma= (diag(1,4)+matrix(1,nrow=4,ncol=4)))
  if(nonlinear){
    y <- rnorm(n, betas[1]*cos(X[,1]-1.5) +betas[2]*X[,2]^3+ X[,3:4]%*%betas[3:4],sd=1)
  }else{
    y <- rnorm(n,X%*%betas,sd=1)
  }
  truedf <- df <- data.frame(y,X)

  ## define LODs for x1 and x2 separately for first half and second half of observations
  LODs <- df[,-1]
  LODs[1:(n/2),1] <- LOD1[1];LODs[(n/2) + (1:(n/2)),1] <- LOD1[2]
  LODs[1:(n/2),2] <- LOD2[1];LODs[(n/2) + (1:(n/2)),2] <- LOD2[2]
  LODs[,3] <- LODs[,4] <- -Inf ## uncensored

  ## censor the values of x1 and x2
  Xcens <- X
  Xcens[Xcens<LODs] <- NA #LODs[Xcens<LODs]
  df[,-1] <- Xcens

  ## prepare LODs for function ## i.e. make LOD NA for non-censored values
  LODs[!is.na(Xcens)] <- NA

  return(list(X=data.frame(X), ## true exposure data
              df=df, ## observed/censored dataset
              LODs=LODs)) ##LODs
}
