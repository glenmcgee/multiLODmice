# ## testing the multiLODmice functions
#
#
# library(crch)
# library(censReg)
# library(mvtnorm)
#
# ## generate data
# set.seed(1)
# n <- 20000
# X <- rmvnorm(n,
#              mean=rep(0,4),
#              sigma= (diag(1,4)+matrix(1,nrow=4,ncol=4)))
# colnames(X) <- paste0("x",1:4)
# y <- rnorm(n,X%*%c(1,1,1,1),sd=1)
# df <- data.frame(y,X)
# g <- lm(y~x1+x2+x3+x4,data=df) ## gold standard
#
# LODs <- df[,-1]
# LODs[1:(n/2),1] <- -2;LODs[(n/2) + (1:(n/2)),1] <- -1
# LODs[1:(n/2),2] <- -1;LODs[(n/2) + (1:(n/2)),2] <- -2
# LODs[,3] <- LODs[,4] <- -Inf ## uncensored
#
# ## censor the values of x1 and x2
# Xcens <- X
# Xcens[Xcens<LODs] <- NA #LODs[Xcens<LODs]
# df[,-1] <- Xcens
#
# ## make LOD for non-censored values -Inf
# LODs[!is.na(Xcens)] <- NA #LODs[Xcens!=LODs] <- NA
#
# ## impute data
# kk <- 20
# df_imp <- multiLODmice(data = df, data.lod = LODs,mi.m = kk)
#
# ## check that imputed values are correctly below the LOD
# for(ii in 1:kk){
#   print(mean(df_imp[[ii]]$x1[is.na(df$x1) & seq(1,n)<(n/2)]<(-2)))
#   print(mean(df_imp[[ii]]$x1[is.na(df$x1) & seq(1,n)>(n/2)]<(-1)))
#   print(mean(df_imp[[ii]]$x2[is.na(df$x2) & seq(1,n)<(n/2)]<(-1)))
#   print(mean(df_imp[[ii]]$x2[is.na(df$x2) & seq(1,n)>(n/2)]<(-2)))
# }
#
# ## compare imputations to true values
# mean_imp1 <- rep(0,sum(is.na(df$x1)))
# mean_imp2 <- rep(0,sum(is.na(df$x2)))
# for(ii in 1:kk){
#   mean_imp1 <- mean_imp1+df_imp[[ii]]$x1[is.na(df$x1)]/kk
#   mean_imp2 <- mean_imp2+df_imp[[ii]]$x2[is.na(df$x2)]/kk
# }
# summary(lm(X[is.na(df$x1),1]~mean_imp1)) ## average imputation looks quite like the true values!
# summary(lm(X[is.na(df$x2),2]~mean_imp2)) ## average imputation looks quite like the true values!
#
# ## Test out MI + Rubin's rules in lm
# ests <- c()
# vcovs <- matrix(0,nrow=5,ncol=5)
# for(ii in 1:kk){
#   g_mi <- lm(y~x1+x2+x3+x4,data=df_imp[[ii]])
#   ests <- rbind(ests,g_mi$coef)
#   vcovs <- vcovs+vcov(g_mi)/kk
# }
# MI_ests <- apply(ests,2,mean) ## works!
# MI_vcov <- vcovs+cov(ests)
#
#
# ############################################################
# ## BY contrast, try with just min LOD (should be invalid)
# minLODs <- LODs
# minLODs$x1[!is.na(minLODs$x1)] <- min(minLODs$x1,na.rm=TRUE)
# minLODs$x2[!is.na(minLODs$x2)] <- min(minLODs$x2,na.rm=TRUE)
#
#
# df_min <- multiLODmice(data = df, data.lod = minLODs,mi.m = kk)
#
# ## compare minutations to true values
# mean_min1 <- rep(0,sum(is.na(df$x1)))
# mean_min2 <- rep(0,sum(is.na(df$x2)))
# for(ii in 1:kk){
#   mean_min1 <- mean_min1+df_min[[ii]]$x1[is.na(df$x1)]/kk
#   mean_min2 <- mean_min2+df_min[[ii]]$x2[is.na(df$x2)]/kk
# }
# summary(lm(X[is.na(df$x1),1]~mean_min1)) ## doesnt look great (as expected)
# summary(lm(X[is.na(df$x2),2]~mean_min2)) ## doesnt look great (as expected)
#
# ## Test out MI + Rubin's rules in lm
# ests <- c()
# vcovs <- matrix(0,nrow=5,ncol=5)
# for(ii in 1:kk){
#   g_mi <- lm(y~x1+x2+x3+x4,data=df_min[[ii]])
#   ests <- rbind(ests,g_mi$coef)
#   vcovs <- vcovs+vcov(g_mi)/kk
# }
# min_ests <- apply(ests,2,mean) ## estimates of beta1,beta2 are biased downwards (as expected)
# min_vcov <- vcovs+cov(ests)
#
#
# ############################################################
# ## BY contrast, try with just max LOD (should be inefficient)
# maxLODs <- LODs
# maxLODs$x1[!is.na(maxLODs$x1)] <- max(maxLODs$x1,na.rm=TRUE)
# maxLODs$x2[!is.na(maxLODs$x2)] <- max(maxLODs$x2,na.rm=TRUE)
#
#
# df_max <- multiLODmice(data = df, data.lod = maxLODs,mi.m = kk)
#
# ## compare maxutations to true values
# mean_max1 <- rep(0,sum(is.na(df$x1)))
# mean_max2 <- rep(0,sum(is.na(df$x2)))
# for(ii in 1:kk){
#   mean_max1 <- mean_max1+df_max[[ii]]$x1[is.na(df$x1)]/kk
#   mean_max2 <- mean_max2+df_max[[ii]]$x2[is.na(df$x2)]/kk
# }
# summary(lm(X[is.na(df$x1),1]~mean_max1)) ## doesnt look great (as expected)
# summary(lm(X[is.na(df$x2),2]~mean_max2)) ## doesnt look great (as expected)
#
# ## Test out MI + Rubin's rules in lm
# ests <- c()
# vcovs <- matrix(0,nrow=5,ncol=5)
# for(ii in 1:kk){
#   g_mi <- lm(y~x1+x2+x3+x4,data=df_max[[ii]])
#   ests <- rbind(ests,g_mi$coef)
#   vcovs <- vcovs+vcov(g_mi)/kk
# }
# maxests <- apply(ests,2,mean) ## estimates of beta1,beta2 are biased downwards (as expected)
# maxvcov <- vcovs+cov(ests)

#
#
# beta <- c(-1,2)
# sig <- 1.5
# n <- 20000
# x <- rnorm(n)
# y <- rnorm(n,beta[1]+beta[2]*x,sig)
# yc <- rcnorm(n,beta[1]+beta[2]*x,sig,left=-2)
# # yt <- y ## by hand gives same thing
# # yt[yt < (-2)] <- -2
# cens <- rep(-2,n)
# df <- data.frame(y,yc,x,cens)
#
#
#
#
#
# #rcnorm(10000,0,1,left=1)
# # gives correct distribution
# ## to see:
# # test <- rcnorm(10000,0,1,left=-0.5)
# # round(mean(test),1) ## only 0 under no truncation
# # round(sd(test),1) ## only 1 under no truncation
# # round(median(test),1) ## same no matter what level of truncation
# # round(quantile(test,0.9),1) ## same no matter what level of truncation
# # hist(test)
#
# ## testing out the truncated regression models
# beta <- c(-1,2)
# sig <- 1.5
# n <- 20000
# x <- rnorm(n)
# y <- rnorm(n,beta[1]+beta[2]*x,sig)
# yc <- rcnorm(n,beta[1]+beta[2]*x,sig,left=-2)
# # yt <- y ## by hand gives same thing
# # yt[yt < (-2)] <- -2
# cens <- rep(-2,n)
# df <- data.frame(y,yc,x,cens)
#
# ## same results under no censoring
# g <- lm(y~x) # same
# g_crch <- crch(y~x,data=df, dist = "gaussian") # same
#
# ## crch works under censoring with single value #lm fails under censoring
# gc <- lm(yc~x) # same
# gc_crch <- crch(yc~x,data=df, dist = "gaussian",left=-2,truncated=FALSE) # same (CENSORED, NOT truncated)
# gc_censReg <- censReg(yc~x,data=df,left=-2)
#
# round(coef(gc)-beta,2)
# round(coef(gc_crch)[1:2]-beta,2)
# round(coef(gc_censReg)[1:2]-beta,2) ## censReg also works with single value
#
# exp(coef(gc_crch)[3]) # gives sigmahat (not sigma2hat)
# exp(coef(gc_censReg)[3]) # censReg also works!
#
#
# ##### MULTIPLE CENSORING VALUES:
# y <- rnorm(n,beta[1]+beta[2]*x,sig)
# yc <- c(rcnorm(n/2,beta[1]+beta[2]*x[1:(n/2)],sig,left=-1),
#         rcnorm(n/2,beta[1]+beta[2]*x[(n/2)+(1:(n/2))],sig,left=-5))
#
# cens <- c(rep(-1,n/2),rep(-5,n/2))
# df <- data.frame(y,yt,yc,x,cens)
#
#
#
# ## crch works under censoring with multiple
# df$ycflat <- df$yc; df$ycflat[df$ycflat< (-1)] <- -1
# gc_crch1 <- crch(ycflat~x,data=df, dist = "gaussian",left=-1,truncated=FALSE) # same (CENSORED, NOT truncated)
# round(coef(gc_crch1)[1:2]-beta,2)
# gc_crch2 <- crch(yc~x,data=df, dist = "gaussian",left=cens,truncated=FALSE) # same (CENSORED, NOT truncated)
# round(coef(gc_crch2)[1:2]-beta,2)
# gc_censReg1 <- censReg(ycflat~x,data=df,left=-1) # max LOD
# round(coef(gc_censReg1)[1:2]-beta,2)
# # gc_censReg2 <- censReg(yc~x,data=df,left=cens) # multiple LODs DOESNT WORK WITH CENSREG
# # round(coef(gc_censReg2)[1:2]-beta,2)
#
# summary(gc_crch1) ## similar estimates
# summary(gc_crch2) # but the one with multiple censoring points has smaller SEs, so suggests it actually is working
#
# ## note that left=-5 would be wrong, since there are other censored values
# gc_crch5 <- crch(yc~x,data=df, dist = "gaussian",left=-5,truncated=FALSE) # same (CENSORED, NOT truncated)
# round(coef(gc_crch5)[1:2]-beta,2)
#
#
# #
# #### Specialized software to do just that!
# library(doMIsaul)
# beta <- c(-1,2)
# sig <- 1.5
# n <- 20000
# thresh <- -2
# x <- rnorm(n)
# y <- rnorm(n,beta[1]+beta[2]*x,sig)
# yc <- rcnorm(n,beta[1]+beta[2]*x,sig,left=thresh)
# df <- data.frame(yc,x)
# Censored <- df[,1,drop=F]
# Censored[Censored$yc!=-2,] <- 1000+runif(sum(Censored$yc!=-2))
# LODs <- Censored[1,,drop=F]
# LODs$yc <- thresh # done
# # LODs$x <- NA
# df2 <- df # done
# df2$yc[df2$yc<=thresh] <- NA # done
# df2$x[10] <- NA # other missingness # done
# test <- MImpute_lcens(data = df2, data.lod = Censored, standards = LODs,
#                       mi.m = 5, mice.log = FALSE)
#
# ## seems to work!
# g1 <- lm(yc~x,data=test[[1]]) ## works
# coef(g1)
# g0 <- lm(yc~x,data=df) ## fails, since it ignores censoring
# coef(g0)
# g2 <- lm(yc~x,data=df2) ## fails, since it is a CC analysis
# coef(g2)
#
# # ### setting same basis matrix for b-splines:
# # g <- gam(tt4 ~ s(triclosan,bs="bs",k=10) + age + bmi + white + urine.sg,data=Xmat,method="REML",family=fam)
# # get_smooth(g,"triclosan")$knots
# # get_smooth(g,"triclosan")$m
# # length(get_smooth(g,"triclosan")$knots ) ## k + m[1]+1
# #
# # FIXEDKNOTS <- list(triclosan=get_smooth(g,"triclosan")$knots )
# # g2 <- gam(tt4 ~ s(triclosan,bs="bs",k=10) + age + bmi + white + urine.sg,
# #           knots=FIXEDKNOTS,
# #           data=Xmat,method="REML",family=fam)
# # get_smooth(g2,"triclosan")$knots
# # get_smooth(g2,"triclosan")$m
#
#
#
#
#
# # #### Specialized software to do just that!
# # # mice.impute.cens will allow the mice function to use the custom imputation method for censored normal data
# # # install.packages("doMIsaul")
# # library(doMIsaul)
#
# ## MImpute_lcens
# toy <- iris[, 1:4]
# # censor on variables 3 and 4, with LOD at quantile .1 and .2.
# LODs <- toy[1, ]
# LODs[1, ] <-c(NA, NA, quantile(toy[,3], .2), quantile(toy[,4], .1))
# # Censor indicator
# Censored <- data.frame(Petal.Length = runif(150, 50,60),
#                        Petal.Width = runif(150, 50,60))
# Censored[toy[,3] < LODs[1, 3], 1] <- LODs[1, 3]
# Censored[toy[,4] < LODs[1, 4], 2] <- LODs[1, 4]
# # NA for censored data
# toy[toy[,3] < LODs[1, 3], 3] <- NA
# toy[toy[,4] < LODs[1, 4], 4] <- NA
# # # Additional missing data
# # toy[sample(1:nrow(toy), 30), 1] <- NA
# # toy[sample(1:nrow(toy), 30), 3] <- NA
# # toy[sample(1:nrow(toy), 30), 4] <- NA
#
# toy.imp <- MImpute_lcens(data = toy, data.lod = Censored, standards = LODs,
#                          mi.m = 5, mice.log = FALSE)
#
#
#
#
#
