#'@title BIC for space under conditional chain graph likelihood setting.
#'@description Calculate the BIC under the conditional chain graph likelihood.
#'@param fit The fitted model object returned from spacemap::[spacemap, space]. 
#'@param NN An integer, the sample size. 
#'@param tol A numeric lower bound on non-zero parameter values.  
#'@param Y A numeric matrix of N X Q responses.
#'@param X A numeric matrix of N X P predictors.
#'@return the Bayesian Information Criterion (BIC) under conditional chain graph likelihood formulation.
#'@references  Zhang L, Kim S (2014) Learning Gene Networks under SNP Perturbations Using eQTL Datasets. 
#'PLoS Comput Biol 10(2): e1003420. doi: 10.1371/journal.pcbi.1003420 
bicCondLike <- function(fit, tol, Y, X, model = c("spacemap", "space", "scggm"), Yindex = NULL, Xindex = NULL) {
  
  if (model == "spacemap") {
    ThetaYY <- calcThetaYY(fit$ParCor, fit$sig.fit)
    ThetaXY <- calcThetaXY(fit$Gamma, fit$sig.fit)
  } else if(!is.null(Yindex) & !is.null(Xindex) & model == "space") {
    stopifnot(is.integer(Yindex), is.integer(Xindex))
    ThetaXXYY <- calcThetaYY(fit$ParCor, fit$sig.fit)
    ThetaYY <- ThetaXXYY[Yindex,Yindex]
    ThetaXY <- ThetaXXYY[Xindex,Yindex]
  } else if (model == "scggm") {
    ThetaYY <- fit$ThetaYY
    ThetaXY <- fit$ThetaXY
  } else {
    stop("bicCondLike: No suitable method selected. Check arguments.")
  }
  
  #conditional gaussian log likelihood
  logLike <- gaussCondLogLik(Y, X, ThetaYY, ThetaXY)
  #degrees of freedom: don't penalize diagonal of YY
  dfParCor <- nonZeroUpper(ThetaYY, tol)
  dfGamma <- nonZeroWhole(ThetaXY, tol) 
  df <- dfParCor + dfGamma
  NN <- nrow(Y)
  bic <- -2*logLike + log(NN)*df
  return(c(bic = bic, logLike = logLike, df = df, dfGamma = dfGamma, dfParCor = dfParCor))
}

bicRSS <- function(model, fit, NN, df = NULL, tol = NULL) {
  if (is.null(df)) {
    if (model == "spacemap") {
      df <- nonZeroWhole(fit$Gamma, tol) + nonZeroUpper(fit$ParCor, tol)
    } else if(model == "space") {
      df <- nonZeroUpper(fit$ParCor, tol)
    } else if (model == "remMap") {
      df <- nonZeroWhole(fit$Gamma, tol)
    } else {
      stop("bicRSS: No suitable method selected for deg. of freedom calculation. Check arguments.")
    }
  }
  logLike <-  NN*log( fit$rss / NN)
  bic <- logLike + df*log(NN)
  return(c(bic = bic, logLike = logLike, df = df))
}

#'@title BIC for space under likelihood setting.
#'@description Calculate the BIC.
#'@param fit The returned object from spacemap::space. 
#'@param NN An integer, the sample size. 
#'@param tol A numeric lower bound on non-zero parameter values.  
#'@return the Bayesian Information Criterion (BIC).
bicLike <- function(fit, tol, dat) {
  ThetaYY <- calcThetaYY(fit$ParCor, fit$sig.fit)
  logLike <- gaussLogLik(ThetaYY, dat)
  df <- nonZeroUpper(fit$ParCor, tol)
  NN <- nrow(dat)
  bic <- NN * (-2 * logLike ) + df*log(NN)
  return(c(bic = bic, logLike = logLike, df = df))
}
