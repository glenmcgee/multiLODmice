## generate complete data
n <- 5000
dat <- gendat_multiLOD(n,nonlinear=TRUE)
X <- dat$X
df <- dat$df
LODs <- dat$LODs

## model fit with no censoring
g_gold <- mgcv::gam(df$y~s(X1)+s(X2)+X3+X4,data=X) ## gold standard
plot(g_gold)


## complete cases only
g_cc <- mgcv::gam(df$y~s(X1)+s(X2)+X3+X4,data=df)

## naive imputation: LOD/sqrt(2)
df_naive <- df
df_naive$X1[!is.na(LODs$X1)] <- LODs$X1[!is.na(LODs$X1)]/sqrt(2)
df_naive$X2[!is.na(LODs$X2)] <- LODs$X2[!is.na(LODs$X2)]/sqrt(2)
g_naive <- mgcv::gam(df$y~s(X1)+s(X2)+X3+X4,data=df_naive)

## MI
kk <- 10
df_imp <- multiLODmice(data = df, data.lod = LODs,mi.m = kk,maxit=30)
g_imp1 <- mgcv::gam(df$y~s(X1)+s(X2)+X3+X4,data=df_imp[[1]],method="REML")
plot(g_imp1)



#### does splines on y help? marginally?
ys <- data.frame(ns(df$y,4));colnames(ys) <- paste0("y",1:4)
df_imp <- multiLODmice(data = data.frame(ys,df[,-1]), data.lod = LODs,mi.m = kk)
g_test <- mgcv::gam(df$y~s(X1)+s(X2)+X3+X4,data=data.frame(y=df$y,df_imp[[1]]),method="REML")
plot(g_test)


#### does LESS info on y help? possibly..but only marginally
df_imp <- multiLODmice(data = data.frame(df[,-1]), data.lod = LODs,mi.m = kk)
g_test <- mgcv::gam(df$y~s(X1)+s(X2)+X3+X4,data=data.frame(y=df$y,df_imp[[1]]),method="REML")
plot(g_test)

#### what about adding other x variables -- QUITE GOOD!!
xs <- data.frame(ns(df$X1,5));colnames(xs) <- paste0("x1s",1:5)
df_imp <- multiLODmice(data = data.frame(xs[,-1],df), data.lod = LODs,mi.m = kk)
g_test <- mgcv::gam(df$y~s(X1)+s(X2)+X3+X4,data=data.frame(y=df$y,df_imp[[1]]),method="REML")
plot(g_test)

#### what about adding both x and y variables --- back to not great
xs <- data.frame(ns(df$X1,5));colnames(xs) <- paste0("x1s",1:5)
ys <- data.frame(ns(df$y,4));colnames(ys) <- paste0("y",1:4)
df_imp <- multiLODmice(data = data.frame(ys,xs[,-1],df[,-1]), data.lod = LODs,mi.m = kk)
g_test <- mgcv::gam(df$y~s(X1)+s(X2)+X3+X4,data=data.frame(y=df$y,df_imp[[1]]),method="REML")
plot(g_test)

#### what about adding other x variables and NO Y ---
xs <- data.frame(ns(df$X1,5));colnames(xs) <- paste0("x1s",1:5)
df_imp <- multiLODmice(data = data.frame(xs[,-1],df[,-1]), data.lod = LODs,mi.m = kk)
g_test <- mgcv::gam(df$y~s(X1)+s(X2)+X3+X4,data=data.frame(y=df$y,df_imp[[1]]),method="REML")
plot(g_test)

plot(g_gold,xlim=c(-4,4),ylim=c(-6,6))
plot(mgcv::gam(df$y~s(X1)+s(X2)+X3+X4,data=data.frame(y=df$y,df_imp[[4]]),method="REML"),xlim=c(-4,4),ylim=c(-6,6))
plot(mgcv::gam(df$y~s(X1)+s(X2)+X3+X4,data=data.frame(y=df$y,df_imp[[9]]),method="REML"),xlim=c(-4,4),ylim=c(-6,6))






### SOME TESTS about the imputation model
n <- 200000
dat <- gendat_multiLOD(n,rho=0,betas=c(2,0.08,1,1),nonlinear=TRUE)
X <- dat$X
df <- dat$df
LODs <- dat$LODs

dftest <- df
dftest[,-1] <- X
dftest$X1gold <- dftest$X1
dftest$X1 <- df$X1

## gold
g_gold <- gam(y~s(X1gold)+s(X2)+X3+X4,data=dftest)
plot(g_gold)


## linear is no good not great
m_imp <- lm(X1~X2+X3+X4+y,data=dftest[!is.na(dftest$X1),])
preds <- predict(m_imp,newdata = dftest[is.na(dftest$X1),])
plot(dftest$X1gold[is.na(dftest$X1)]~preds)
# plot(dftest$X1gold[!is.na(dftest$X1)]~y[!is.na(dftest$X1)])


#### DEPENDING ON Y
## m observed
m_obs <- gam(X1~X2+X3+X4+s(y),data=dftest[!is.na(dftest$X1),])
pred_obs <- predict(m_obs,newdata = dftest[is.na(dftest$X1),]) ## predict for unobserved
summary(lm(dftest$X1gold[is.na(dftest$X1)]~pred_obs))
plot(m_obs)

## unobserved
m_unobs <- gam(X1gold~X2+X3+X4+s(y),data=dftest[is.na(dftest$X1),])
pred_unobs <- predict(m_unobs,newdata = dftest[is.na(dftest$X1),]) ## here these are basically fitted vals
summary(lm(dftest$X1gold[is.na(dftest$X1)]~pred_unobs))
plot(m_unobs)

## both
m_both <- gam(X1gold~X2+X3+X4+s(y),data=dftest)
plot(m_both)
pred_both <- predict(m_both,newdata = dftest[is.na(dftest$X1),]) ## predict for unobserved
summary(lm(dftest$X1gold[is.na(dftest$X1)]~pred_both))
plot(dftest$X1gold[is.na(dftest$X1)]~pred_both)

### THE RELATIONSHIP AMONG OBSERVED IS FUNDAMENTALLY DIFFERENT THAN AMONG UNOBSERVED


#### NO Y
## m observed
m_obs <- gam(X1~X2+X3+X4,data=dftest[!is.na(dftest$X1),])
pred_obs <- predict(m_obs,newdata = dftest[is.na(dftest$X1),]) ## predict for unobserved
summary(lm(dftest$X1gold[is.na(dftest$X1)]~pred_obs))
plot(dftest$X1gold[is.na(dftest$X1)]~pred_obs)

## unobserved
m_unobs <- gam(X1gold~X2+X3+X4,data=dftest[is.na(dftest$X1),])
pred_unobs <- predict(m_unobs,newdata = dftest[is.na(dftest$X1),]) ## here these are basically fitted vals
summary(lm(dftest$X1gold[is.na(dftest$X1)]~pred_unobs))
plot(dftest$X1gold[is.na(dftest$X1)]~pred_unobs)

## both
m_both <- gam(X1gold~X2+X3+X4,data=dftest)
pred_both <- predict(m_both,newdata = dftest[is.na(dftest$X1),]) ## predict for unobserved
summary(lm(dftest$X1gold[is.na(dftest$X1)]~pred_both))
plot(dftest$X1gold[is.na(dftest$X1)]~pred_both)

### OF COURSE THEY ALL PERFORM THE SAME--WE'VE SIMULATED UNDER INDEPENDENT X's
### JUST BY CHANGE WE DO BETTER THAN ABOVE WHERE WE DEFINITELY ASSUME THE WRONG FUNCTIONAL RELATIONSHIP



### WHAT IF NOW WE HAVE DEPENDENCE AMONG X's???

n <- 200000
dat <- gendat_multiLOD(n,rho=0.7,betas=c(2,0.08,1,1),nonlinear=TRUE)
X <- dat$X
df <- dat$df
LODs <- dat$LODs

dftest <- df
dftest[,-1] <- X
dftest$X1gold <- dftest$X1
dftest$X1 <- df$X1


#### DEPENDING ON Y
## m observed
m_obs <- gam(X1~X2+X3+X4+s(y),data=dftest[!is.na(dftest$X1),])
pred_obs <- predict(m_obs,newdata = dftest[is.na(dftest$X1),]) ## predict for unobserved
summary(lm(dftest$X1gold[is.na(dftest$X1)]~pred_obs))
plot(m_obs)

## unobserved
m_unobs <- gam(X1gold~X2+X3+X4+s(y),data=dftest[is.na(dftest$X1),])
pred_unobs <- predict(m_unobs,newdata = dftest[is.na(dftest$X1),]) ## here these are basically fitted vals
summary(lm(dftest$X1gold[is.na(dftest$X1)]~pred_unobs))
plot(m_unobs)

## both
m_both <- gam(X1gold~X2+X3+X4+s(y),data=dftest)
plot(m_both)
pred_both <- predict(m_both,newdata = dftest[is.na(dftest$X1),]) ## predict for unobserved
summary(lm(dftest$X1gold[is.na(dftest$X1)]~pred_both))
plot(dftest$X1gold[is.na(dftest$X1)]~pred_both)

### THE RELATIONSHIP AMONG OBSERVED IS FUNDAMENTALLY DIFFERENT THAN AMONG UNOBSERVED


#### NO Y
## m observed
m_obs <- gam(X1~X2+X3+X4,data=dftest[!is.na(dftest$X1),])
pred_obs <- predict(m_obs,newdata = dftest[is.na(dftest$X1),]) ## predict for unobserved
summary(lm(dftest$X1gold[is.na(dftest$X1)]~pred_obs))
plot(dftest$X1gold[is.na(dftest$X1)]~pred_obs)

## unobserved
m_unobs <- gam(X1gold~X2+X3+X4,data=dftest[is.na(dftest$X1),])
pred_unobs <- predict(m_unobs,newdata = dftest[is.na(dftest$X1),]) ## here these are basically fitted vals
summary(lm(dftest$X1gold[is.na(dftest$X1)]~pred_unobs))
plot(dftest$X1gold[is.na(dftest$X1)]~pred_unobs)

## both
m_both <- gam(X1gold~X2+X3+X4,data=dftest)
pred_both <- predict(m_both,newdata = dftest[is.na(dftest$X1),]) ## predict for unobserved
summary(lm(dftest$X1gold[is.na(dftest$X1)]~pred_both))
plot(dftest$X1gold[is.na(dftest$X1)]~pred_both)

### OF COURSE THEY ALL PERFORM THE SAME--WE'VE SIMULATED UNDER INDEPENDENT X's
### JUST BY CHANGE WE DO BETTER THAN ABOVE WHERE WE DEFINITELY ASSUME THE WRONG FUNCTIONAL RELATIONSHIP


