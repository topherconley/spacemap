

library(spacemap)
data(sim1)
dat <- sim1[c("XY", "X", "Y", "Xindex", "Yindex")]
N <- nrow(dat$X)
P <- ncol(dat$X)
Q <- ncol(dat$Y)
tol <- 1e-6

trueParCor <- sim1$trueParCor


#spacemap

#equal non-zero entries between fit and refit
fit <- spacemap(Y.m = dat$Y, X.m = dat$X, slasso = 69.67802 , rlasso = 28.84801 , rgroup = 12.38648)
nzfitpc <- (abs(fit$ParCor) > tol) + 0
nzfitg <- (abs(fit$Gamma) > tol) + 0
sparse <- olsRefitSpacemap(Y = dat$Y, X = dat$XY, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma = fit$sig.fit, RSS = fit$rss, tol = tol, ridge = 0)
nzrefitpc <- (abs(fit$ParCor) > tol) + 0
nzrefitg <- (abs(fit$Gamma) > tol) + 0
stopifnot(identical(nzfitpc, nzrefitpc), identical(nzfitg, nzrefitg))


#How different is fastLmPure and my implementation of refit
fit <- spacemap(Y.m = dat$Y, X.m = dat$X, slasso = 69.67802 , rlasso = 28.84801 , rgroup = 12.38648)
nzfit_pc_val <- fit$ParCor[as.logical(nzfitpc)]
nzfit_g_val <- fit$Gamma[as.logical(nzfitg)]

#r refit
rrefit <- OLS.Sol2(Y = dat$Y, dat$X, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma.fit = fit$sig.fit, tol = tol)
nzrefit1_pc_val <- rrefit$ParCor[as.logical(nzfitpc)]
nzrefit1_g_val <- rrefit$Gamma[as.logical(nzfitg)]
sparse <- olsRefitSpacemap(Y = dat$Y, X = dat$XY, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma = fit$sig.fit, RSS = fit$rss, tol = tol, ridge = 0)
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
fit <- spacemap(Y.m = dat$Y, X.m = dat$X, slasso = 69.67802 , rlasso = 28.84801 , rgroup = 12.38648)
fit$rss
all.equal(trueParCor$yy, fit$ParCor)
all.equal(trueParCor$xy, fit$Gamma)

rrefit <- OLS.Sol2(Y = dat$Y, dat$X, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma.fit = fit$sig.fit, tol = tol)
olsRSSSpacemap(Y = dat$Y, X = dat$X, refit = rrefit)
all.equal(trueParCor$yy, rrefit$ParCor)
all.equal(trueParCor$xy, rrefit$Gamma)




library(microbenchmark)
microbenchmark(OLS.Sol2(Y = dat$Y, dat$X, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma.fit = fit$sig.fit, tol = tol), 
               olsRefitSpacemap(Y = dat$Y, X = dat$XY, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma = fit$sig.fit, RSS = fit$rss, tol = tol, ridge = 0))


refit <- olsRefitSpacemap(Y = dat$Y, X = dat$XY, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma = fit$sig.fit, RSS = fit$rss, tol = tol, ridge = 0)
fit$rss
all.equal(trueParCor$yy, fit$ParCor)
all.equal(trueParCor$xy, fit$Gamma)


all.equal(rrefit$ParCor, fit$ParCor)






#refit has lower rss
olsRSSSpacemap <- function(Y, X, refit) { 
  diag(refit$ParCor) <- 0
  Y.pred <- Y%*%refit$ParCor + X%*%refit$Gamma
  resid <- Y  - Y.pred
  rss <- sum(resid*resid) 
  rss
}

fit <- spacemap(Y.m = dat$Y, X.m = dat$X, slasso = 69.67802 , rlasso = 28.84801 , rgroup = 12.38648)
fit$rss
rrefit <- OLS.Sol2(Y = dat$Y, dat$X, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma.fit = fit$sig.fit, tol = tol)
olsRSSSpacemap(Y = dat$Y, X = dat$X, refit = rrefit)
refit <- olsRefitSpacemap(Y = dat$Y, X = dat$XY, ParCor = fit$ParCor, Gamma = fit$Gamma, sigma = fit$sig.fit, RSS = fit$rss, tol = tol, ridge = 0)
fit$rss
all.equal(rrefit$ParCor, fit$ParCor)

all.equal(trueParCor$xy, fit$Gamma)
all.equal(trueParCor$xy, rrefit$Gamma)


all.equal(trueParCor$yy, fit$ParCor)
all.equal(trueParCor$yy, rrefit$ParCor)

