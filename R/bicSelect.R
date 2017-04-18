# 
# 
# 
# #'@title BIC for spacemap
# #'@description Calculate the pseudo-likelihood BIC
# #'@param spmap The returned object from spacemap::spacemap. 
# #'@param NN An integer, the sample size. 
# #'@param tol A numeric lower bound on non-zero parameter values.  
# #'@return the Bayesian Information Criterion (BIC) under pseudo-likelihood formulation.
# bicSpaceMap <- function(spmap, NN, tol) {
#   #multiply degree freedom of partial correlations of YY by two 
#   #b/c of BIC-type is thought of as a sum of Q reqressions. 
#   #See space paper for details. 
#   degfree <- sum(abs(spmap$Gamma) > tol) + 2*sum(abs(spmap$ParCor[upper.tri(spmap$ParCor)]) > tol)
#   NN*log( spmap$rss / NN) + degfree*log(NN)
# }
# 
# #'@title BIC for remMap
# #'@description Calculate the pseudo-likelihood BIC
# #'@param rmap The returned object from spacemap::remMap. 
# #'@param NN An integer, the sample size. 
# #'@param tol A numeric lower bound on non-zero parameter values.  
# #'@return the Bayesian Information Criterion (BIC) under pseudo-likelihood formulation.
# bicRemMap <- function(rmap, NN, tol) {
#   degfree <- sum(abs(rmap$Gamma) > tol) 
#   NN*log( rmap$rss / NN) + degfree*log(NN)
# }
# 
# #'@title BIC for space
# #'@description Calculate the pseudo-likelihood BIC
# #'@param sp The returned object from spacemap::space. 
# #'@param NN An integer, the sample size. 
# #'@param tol A numeric lower bound on non-zero parameter values.  
# #'@return the Bayesian Information Criterion (BIC) under pseudo-likelihood formulation.
# #'@
# bicSpace <- function(sp, NN, tol) {
#   #multiply degree freedom of partial correlations of YY by two 
#   #b/c of BIC-type is thought of as a sum of Q reqressions. 
#   #See space paper for details. 
#   degfree <- 2*sum(abs(sp$ParCor[upper.tri(sp$ParCor)]) > tol)
#   NN*log( sp$rss / NN) + degfree*log(NN)
# }
# 
# #'@title BIC for space under likelihood setting.
# #'@description Calculate the pseudo-likelihood BIC.
# #'@param sp The returned object from spacemap::space. 
# #'@param NN An integer, the sample size. 
# #'@param tol A numeric lower bound on non-zero parameter values.  
# #'@return the Bayesian Information Criterion (BIC) under pseudo-likelihood formulation.
# bicLikeSpace <- function(sp, NN, tol, dat) {
#   degfree <- sum(abs(sp$ParCor)>1e-6)
#   p <-ncol(sp$ParCor)
#   sqrt_dii_djj <- matrix(sqrt(sp$sig.fit),p,p,byrow=FALSE)*matrix(sqrt(sp$sig.fit),p,p,byrow=TRUE)
#   Sigma.fit.inv <- -1 * sp$ParCor  *  sqrt_dii_djj
#   diag(Sigma.fit.inv) <- sp$sig.fit
#   #Sample covariance matrix
#   S <- cov(dat)
#   #BIC 
#   NN*( -1*log(det(Sigma.fit.inv)) + sum(diag(Sigma.fit.inv%*%S)) ) + degfree*log(NN)
# }
# 
# #'@title BIC for space under conditional chain graph likelihood setting.
# #'@description Calculate the BIC under the conditional chain graph likelihood.
# #'@param fit The fitted model object returned from spacemap::[spacemap, space]. 
# #'@param NN An integer, the sample size. 
# #'@param tol A numeric lower bound on non-zero parameter values.  
# #'@param Y A numeric matrix of N X Q responses.
# #'@param X A numeric matrix of N X P predictors.
# #'@return the Bayesian Information Criterion (BIC) under conditional chain graph likelihood formulation.
# #'@references  Zhang L, Kim S (2014) Learning Gene Networks under SNP Perturbations Using eQTL Datasets. 
# #'PLoS Comput Biol 10(2): e1003420. doi: 10.1371/journal.pcbi.1003420 
# bicCondLike <- function(fit, NN, tol, Y, X, model = c("spacemap", "space", "scggm"), Yindex = NULL, Xindex = NULL) {
#  
#   if (model == "spacemap") {
#     ThetaYY <- calcThetaYY(fit$ParCor, fit$sig.fit)
#     ThetaXY <- calcThetaXY(fit$Gamma, fit$sig.fit)
#   } else if(!is.null(Yindex) & !is.null(Xindex) & model == "space") {
#     stopifnot(is.integer(Yindex), is.integer(Xindex))
#     ThetaXXYY <- calcThetaYY(fit$ParCor, fit$sig.fit)
#     ThetaYY <- ThetaXXYY[Yindex,Yindex]
#     ThetaXY <- ThetaXXYY[Xindex,Yindex]
#   } else if (model == "scggm") {
#     ThetaYY <- fit$ThetaYY
#     ThetaXY <- fit$ThetaXY
#   } else {
#     stop("bicCondLike: No suitable method selected. Check arguments.")
#   }
# 
#   #call RcppArmadillo
#   logLik <- gaussCondLogLik(Y, X, ThetaYY, ThetaXY)
#   
#   #degrees of freedom: don't penalize diagonal of YY
#   degFreeLambda <- sum(abs(ThetaXY) > tol) + sum(abs(ThetaYY[upper.tri(ThetaYY)]) > tol)
#   
#   return(-2*logLik + log(NN)*degFreeLambda)
# }
# # 
# # bicSelect <- function(Y.m, X.m, tune, method = c("spacemap", "space", "remmap"), corenum = NULL, sig = NULL, weight = NULL, remWeight = NULL, 
# #                       iter = 1, tol = 1e-6, cd_iter = 1e7) {
# #   
# #   library(foreach)
# #   library(doParallel)
# #   #parallel functionality
# #   if (!is.null(corenum)) {
# #     cl <- makeCluster(corenum)
# #     registerDoParallel(cl)
# #   }
# #   
# #   NN <- nrow(Y.m)
# #   
# #   bic <- foreach(l = 1:nrow(tune)) %dopar% {
# #     
# #     if (method == "spacemap") {
# #       fit <- spacemap(Y.m = Y.m, X.m = X.m, 
# #                       lamS1 = tune$s1[l], lamS2=0, lamR1 = tune$r1[l],  lamR2 = tune$r2[l], 
# #                       sig = sig, weight = weight, remWeight = remWeight, iter=iter,
# #                       tol = tol, cd_iter = cd_iter)
# #       bic <- bicSpaceMap(fit, NN, tol)
# #     } else if (method == "remmap") {
# #       #run remMap
# #       fit <- remMap(X.m, Y.m, lamL1 = lseq$r1[l], lamL2 = lseq$r2[l], 
# #                     phi0 = NULL, W.m = NULL, tol = tol)
# #       bic <- bicRemMap(fit, NN, tol)
# #     } else if (method == "space") {
# #       fit <- space(cbind(X.m, Y.m), lam1 = lseq$s1[l], lam2=0, iter=iter, 
# #                          cd_iter = cd_iter, tol = tol)
# #       bicSpace(fit, NN, tol)
# #     } else {
# #       bic <- NA
# #       fit <- NA
# #     }
# #     list(bic = bic, fit = fit)
# #   }
# #   
# #   #clean up parallel forks
# #   if (!is.null(corenum))
# #     stopCluster(cl)
# #   
# #   return(list(bic = bic, tune = tune))
# # }
