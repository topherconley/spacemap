
#Fit the space over 

#Arguments
#model type
# tuning grid
# additional arguments

initFit <- function(data, method = c("spacemap", "space"), tuneGrid, refit = TRUE, ...) {
  #additional arguments
  library(foreach)
  opt <- list(...)
  NN <- nrow(data$Y)
  if (method == "spacemap") {
    
    outfits <- foreach(l = seq_len(nrow(tuneGrid)), .combine = 'rbind') %dopar% {
      fit <- spacemap(Y.m = data$Y, X.m = data$X, 
                      slasso = tuneGrid$slasso[l], sridge = opt$sridge, 
                      rlasso = tuneGrid$rlasso[l], rgroup = tuneGrid$rgroup[l], 
                      sig=opt$sig, weight=opt$weight, remWeight = opt$remWeight, iter= opt$iter,
                      tol = opt$tol, cd_iter = opt$cd_iter)
      if (fit$convergence) {
        conv1 <- 1
        
        if (refit) {
          #recompute ols estimates: note change fit object by reference
          sparseFit <- olsRefitSpacemap(Y = data$Y, X = data$X, 
                                        ParCor = fit$ParCor, Gamma = fit$Gamma, 
                                        sigma = fit$sig.fit, RSS = fit$rss, tol = opt$tol)  
        } else {
          sparseFit <- TRUE
        }
    
        
        
        if(sparseFit) {
          conv2 <- 1
          bcl <- bicCondLike(fit = fit, tol = opt$tol, model = method,
                             Y = data$Y, X = data$X)
        } else {
          conv2 <- 0
          bcl <- NA
          brss <- NA
        } # if else refit convergence
      } else {
        conv1 <- 0
        conv2 <- 0
        bcl <- NA
        brss <- NA
      } # if else convergence
      
      #degrees of freedom
      dfParCor <- nonZeroUpper(fit$ParCor, opt$tol)
      dfGamma <- nonZeroWhole(fit$Gamma, opt$tol) 
      degFree <- dfParCor + dfGamma
      brss <- NN*log( fit$rss / NN) + degFree*log(NN)
      
      
      matrix(c(conv1, conv2, bcl["bic"], brss, degFree, dfParCor, dfGamma), nrow = 1, ncol = 7)
    } # end foreach over tuneGrid
    outfits <- as.data.frame(outfits)
    names(outfits) <- c("convFit", "convRefit", "bicCondLike", "bicRSS", "dfTotal", "dfParCor", "dfGamma")
  } else if (method == "space") {
    
    outfits <- foreach(l = seq_len(nrow(tuneGrid)), .combine = 'rbind') %dopar% {
      fit <- spacemap::space.joint(Y.m = data$XY, lam1 = tuneGrid$slasso[l], lam2 = opt$sridge, iter = opt$iter, 
                                   cd_iter = opt$cd_iter, tol = opt$tol)
      if (fit$convergence) {
        conv1 <- 1
        
        if (refit) {
          sparseFit <- spacemap::olsRefitSpace(Y = data$XY,
                                               ParCor = fit$ParCor, 
                                               sigma = fit$sig.fit, RSS = fit$rss, tol = opt$tol) 
          #sparseFit <- TRUE
          #fit <- spacemap::rolsRefit(fit, data$XY, tol = opt$tol)  
        } else { 
          sparseFit <- TRUE
        }
        
        if(sparseFit) {
          conv2 <- 1
          #bl <- bicLike(fit = fit, tol = opt$tol, dat = data$XY)
          #brss <- bicRSS(model = method, fit = fit, NN = NN, df = bl["df"], tol = opt$tol)
          bcl <- bicCondLike(fit = fit, tol = opt$tol, model = method,
                             Y = data$Y, X = data$X,
                             Yindex = opt$Yindex, Xindex = opt$Xindex)
          degFree <- bcl["df"]
          dfGamma <- bcl["dfGamma"]
          dfParCor <- bcl["dfParCor"]
          brss <- NN*log( fit$rss / NN) + degFree*log(NN)
        } else {
          conv2 <- 0
          bcl <- NA
          brss <- NA
          dfGamma <- NA
          dfParCor <-NA
          degFree <- nonZeroUpper(fit$ParCor, opt$tol)
          #        bcl <- NA
        } # if else refit convergence
      } else {
        conv1 <- 0
        conv2 <- 0
        bl <- NA
        brss <- NA
        dfGamma <- NA
        dfParCor <-NA
        degFree <- nonZeroUpper(fit$ParCor, opt$tol)
      }  # if else convergence
      
      matrix(c(conv1, conv2, bcl["bic"], brss, degFree, dfGamma, dfParCor), nrow = 1, ncol = 7)
    }  # end foreach over tuneGrid
    outfits <- as.data.frame(outfits)
    names(outfits) <- c("convFit", "convRefit", "bicCondLike", "bicRSS", "dfTotal", "dfGamma", "dfParCor")
  } else {
    message("Incorrect method arguments.")
    return(NULL)
  }
  outfits
}

#
