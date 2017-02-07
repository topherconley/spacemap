
library(spacemap)
datdir <- "~/scratch-data/sim-spacemap/hub11-20/2016/datasets/n250"
dat <- readRDS(file = file.path(datdir, "hub11-20-250-dataid_002.rds"))
trueParCor <- dat$trueParCor

N <- nrow(dat$X)
P <- ncol(dat$X)
Q <- ncol(dat$Y)
tol <- 1e-6

slasso <- 100
rlasso <- 44
rgroup <- 20

fit <- spacemap(Y.m = dat$Y, X.m = dat$X, slasso = slasso, rlasso = rlasso, rgroup = rgroup)
nonZeroUpper(fit$ParCor, tol = 1e-6)
nonZeroWhole(fit$Gamma, tol = 1e-6)

#spacemap
#equal non-zero entries between fit and refit
fit <- spacemap(Y.m = dat$Y, X.m = dat$X, slasso = slasso, rlasso = rlasso, rgroup = rgroup)
nzfitpc <- (abs(fit$ParCor) > tol) + 0
nzfitg <- (abs(fit$Gamma) > tol) + 0
sparse <- olsRefitSpacemap(Y = dat$Y, X = dat$X, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma = fit$sig.fit, RSS = fit$rss, tol = tol, ridge = 0)
nzrefitpc <- (abs(fit$ParCor) > tol) + 0
nzrefitg <- (abs(fit$Gamma) > tol) + 0
stopifnot(identical(nzfitpc, nzrefitpc), identical(nzfitg, nzrefitg))



#How different is fastLmPure and my implementation of refit
fit <- spacemap(Y.m = dat$Y, X.m = dat$X, slasso = slasso , rlasso = rlasso, rgroup = rgroup)
nzfit_pc_val <- fit$ParCor[as.logical(nzfitpc)]
nzfit_g_val <- fit$Gamma[as.logical(nzfitg)]

#r refit
rrefit <- OLS.Sol2(Y = dat$Y, dat$X, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma.fit = fit$sig.fit, tol = tol)
nzrefit1_pc_val <- rrefit$ParCor[as.logical(nzfitpc)]
nzrefit1_g_val <- rrefit$Gamma[as.logical(nzfitg)]
sparse <- olsRefitSpacemap(Y = dat$Y, X = dat$X, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma = fit$sig.fit, RSS = fit$rss, tol = tol, ridge = 0)
nzrefit2_pc_val <- fit$ParCor[as.logical(nzrefitpc)]
nzrefit2_g_val <- fit$Gamma[as.logical(nzrefitg)]

plot(nzrefit1_pc_val, nzrefit2_pc_val)
abline(0,1)
plot(nzrefit1_g_val, nzrefit2_g_val)
abline(0,1)

plot(nzfit_g_val, nzrefit2_g_val)
abline(0,1)
plot(nzfit_g_val, nzrefit1_g_val)
abline(0,1)

all.equal(nzrefit1_pc_val, nzrefit2_pc_val, tolerance  = tol)
all.equal(nzrefit1_g_val, nzrefit2_g_val, tolerance  = tol)


#less biased implementation
fit <- spacemap(Y.m = dat$Y, X.m = dat$X, slasso = slasso , rlasso = rlasso, rgroup = rgroup)
fit$rss
all.equal(trueParCor$yy, fit$ParCor)
all.equal(trueParCor$xy, fit$Gamma)

rrefit <- OLS.Sol2(Y = dat$Y, dat$X, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma.fit = fit$sig.fit, tol = tol)

#refit has lower rss
olsRSSSpacemap <- function(Y, X, refit) { 
  diag(refit$ParCor) <- 0
  Y.pred <- Y%*%refit$ParCor + X%*%refit$Gamma
  resid <- Y  - Y.pred
  rss <- sum(resid*resid) 
  rss
}

olsRSSSpacemap(Y = dat$Y, X = dat$X, refit = rrefit)
all.equal(trueParCor$yy, rrefit$ParCor)
all.equal(trueParCor$xy, rrefit$Gamma)


refit <- olsRefitSpacemap(Y = dat$Y, X = dat$X, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma = fit$sig.fit, RSS = fit$rss, tol = tol, ridge = 0)
fit$rss

plot(trueParCor$yy, fit$ParCor)
all.equal(trueParCor$yy, fit$ParCor)
plot(trueParCor$xy, fit$Gamma)
abline(0,1)
all.equal(trueParCor$xy, fit$Gamma)
mean((trueParCor$xy - fit$Gamma)^2)

#remMap fit after 
W.m <- matrix(Inf, P, Q)
W.m[abs(fit$Gamma) > tol] <- 0.0
rm <- remMap(X.m = dat$X, Y.m = dat$Y, lamL1 = rlasso, lamL2 = rgroup)
all.equal(trueParCor$xy, rm$Gamma)

#timing
library(microbenchmark)
microbenchmark(OLS.Sol2(Y = dat$Y, dat$X, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma.fit = fit$sig.fit, tol = tol), 
               olsRefitSpacemap(Y = dat$Y, X = dat$X, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma = fit$sig.fit, RSS = fit$rss, tol = tol, ridge = 0))




##########################SPACE################################

Xindex <- dat$Xindex
Yindex <- dat$Yindex
#equal non-zero entries between fit and refit
fit <- space.joint(Y.m = dat$XY, lam1 = slasso, lam2 = 0)
nzfitpc <- (abs(fit$ParCor[Yindex,Yindex]) > tol) + 0
nzfitg <- (abs(fit$ParCor[Xindex,Yindex]) > tol) + 0
sparse <- olsRefitSpace(Y = dat$XY, fit$ParCor, sigma = fit$sig.fit, RSS = fit$rss, tol = tol, ridge = 0)
nzrefitpc <- (abs(fit$ParCor[Yindex,Yindex]) > tol) + 0
nzrefitg <- (abs(fit$ParCor[Xindex,Yindex]) > tol) + 0
stopifnot(identical(nzfitpc, nzrefitpc), identical(nzfitg, nzrefitg))

#How different is fastLmPure and my implementation of refit
fit <- space.joint(Y.m = dat$XY, lam1 = slasso, lam2 = 0)
nzfit_pc_val <- fit$ParCor[Yindex,Yindex]
nzfit_g_val <- fit$ParCor[Xindex,Yindex]

#r refit
# rrefit <- OLS.Sol2(Y = dat$Y, dat$X, ParCor = fit$ParCor[Yindex,Yindex], Gamma = fit$ParCor[Xindex,Yindex], sigma.fit = fit$sig.fit, tol = tol)
# nzrefit1_pc_val <- rrefit$ParCor[as.logical(nzfitpc)]
# nzrefit1_g_val <- rrefit$Gamma[as.logical(nzfitg)]
rrefit <- rolsRefit(spaceObj = fit, XY = dat$XY, tol = tol)
nzrefit1_pc_val <- rrefit$ParCor[Yindex,Yindex]
nzrefit1_g_val <- rrefit$ParCor[Xindex,Yindex]

sparse <- olsRefitSpace(Y = dat$XY, fit$ParCor, sigma = fit$sig.fit, RSS = fit$rss, tol = tol, ridge = 0)
nzrefit2_pc_val <- fit$ParCor[Yindex,Yindex]
nzrefit2_g_val <- fit$ParCor[Xindex,Yindex]

plot(nzrefit1_pc_val, nzrefit2_pc_val)
abline(0,1)
plot(nzrefit1_g_val, nzrefit2_g_val)
abline(0,1)

plot(nzfit_g_val, nzrefit2_g_val)
abline(0,1)
plot(nzfit_g_val, nzrefit1_g_val)
abline(0,1)

all.equal(nzrefit1_pc_val, nzrefit2_pc_val, tolerance  = tol)
all.equal(nzrefit1_g_val, nzrefit2_g_val, tolerance  = tol)


#less biased implementation
#just y
fit <- space.joint(Y.m = dat$Y, lam1 = slasso, lam2 = 0)
all.equal(trueParCor$yy, fit$ParCor)

#y and x
fit <- space.joint(Y.m = dat$XY, lam1 = slasso, lam2 = 0)
fit$rss
all.equal(trueParCor$yy, fit$ParCor[Yindex,Yindex])
all.equal(trueParCor$xy, fit$ParCor[Xindex,Yindex])



# rrefit <- OLS.Sol2(Y = dat$Y, dat$X, ParCor = fit$ParCor[Yindex,Yindex], Gamma = fit$ParCor[Xindex,Yindex], sigma.fit = fit$sig.fit, tol = tol)
# 
# #refit has lower rss
# olsRSSSpacemap <- function(Y, X, refit) { 
#   diag(refit$ParCor) <- 0
#   Y.pred <- Y%*%refit$ParCor + X%*%refit$Gamma
#   resid <- Y  - Y.pred
#   rss <- sum(resid*resid) 
#   rss
# }
#olsRSSSpacemap(Y = dat$Y, X = dat$X, refit = rrefit)


rrefit <- rolsRefit(spaceObj = fit, XY = dat$XY, tol = tol)
rrefit$rss
all.equal(trueParCor$yy, rrefit$ParCor[Yindex,Yindex])
all.equal(trueParCor$xy, rrefit$ParCor[Xindex,Yindex])

sparse <- olsRefitSpace(Y = dat$XY, fit$ParCor, sigma = fit$sig.fit, RSS = fit$rss, tol = tol, ridge = 0)
fit$rss
all.equal(trueParCor$yy, fit$ParCor[Yindex,Yindex])
all.equal(trueParCor$xy, fit$ParCor[Xindex,Yindex])




