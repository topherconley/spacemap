#######################
testSetPartition <- function(n, fold) {
  randOrder <- sample.int(n)
  foldIds <- list()
  foldSize=floor(n/fold)
  foldBase <- 1:foldSize
  for (f in 1:(fold-1)) {  
    foldIds[[f]]<-randOrder[(f-1)*foldSize+foldBase] 
  } 
  foldIds[[fold]]<-randOrder[((fold-1)*foldSize+1):n]
  foldIds
}
#######################
trainModel <- function(train, method, tuneGrid, refit = TRUE, fold_id, tune_id, ...) {
  
  opt <- list(...)
  #hedge against multicolinearity, for numerical stability
  if (is.null(opt$refitRidge)) opt$refitRidge <- 0.001
  
  #feedback on ols refit assumption of model sparsity
  sparseFit <- TRUE
  library(Matrix)
  if (method == "spacemap") {
    fit <- spacemap(Y.m = train$Y, X.m = train$X, 
                    slasso = tuneGrid$slasso, sridge = opt$sridge, 
                    rlasso = tuneGrid$rlasso, rgroup = tuneGrid$rgroup, 
                    sig=opt$sig, weight=opt$weight, remWeight = opt$remWeight, iter= opt$iter,
                    tol = opt$tol, cd_iter = opt$cd_iter, iscale = opt$iscale)
    #zero out those below tolerance. 
    fit$ParCor[abs(fit$ParCor) <= opt$tol] <- 0.0
    fit$Gamma[abs(fit$Gamma) <= opt$tol] <- 0.0
    #reject models that did not converge
    if(!fit$convergence) return(list(fit = fit, sparseFit = FALSE))
    #recompute ols estimates based on training set: note change fit object by reference
    if(refit) {
      sparseFit <- olsRefitSpacemap(Y = train$Y, X = train$X, 
                                    ParCor = fit$ParCor, Gamma = fit$Gamma, 
                                    sigma = fit$sig.fit, RSS = fit$rss, tol = opt$tol, 
                                    ridge = opt$refitRidge)
    }
    if (!is.null(opt$res_path)) { 
      out <- list(ThetaXY = Matrix(fit$Gamma), ThetaYY = Matrix(fit$ParCor), sig = fit$sig.fit,
                  tune_id = tune_id, fold_id = fold_id, convergence = fit$convergence, sparseFit = sparseFit)
      #"/home/cconley/scratch-data/cptac/tuning/cvfits/spacemap_cptac_rna_cv_tuneid_"
      outfile <- file.path(opt$res_path, paste0("tuneid_", sprintf("%03d", tune_id), "_fold_", sprintf("%02d", fold_id), ".rds"))
      saveRDS(out, file = outfile)
    }
  } else if (method == "space") {
    fit <- space.joint(Y.m = train$XY, lam1 = tuneGrid, lam2 = opt$sridge, iter = opt$iter, 
                       cd_iter = opt$cd_iter, tol = opt$tol, iscale = opt$iscale,
                       sig = opt$sig.fit, rho = opt$ParCor)
    #zero out those below tolerance. 
    fit$ParCor[abs(fit$ParCor) <= opt$tol] <- 0.0
    if(!fit$convergence) return(list(fit = fit, sparseFit = FALSE))
    if(refit) {
      sparseFit <- olsRefitSpace(Y = train$XY,
                                 ParCor = fit$ParCor, 
                                 sigma = fit$sig.fit, RSS = fit$rss, tol = opt$tol,
                                 ridge = opt$refitRidge)
    }
    if (!is.null(opt$res_path)) { 
      out <- list(ThetaXY = Matrix(fit$ParCor[opt$Xindex,opt$Yindex]), 
                  ThetaYY = Matrix(fit$ParCor[opt$Yindex,opt$Yindex]),
                  sig = fit$sig.fit,
                  tune_id = tune_id, fold_id = fold_id, convergence = fit$convergence, sparseFit = sparseFit)
      #"/home/cconley/scratch-data/cptac/tuning/cvfits/space_cptac_rna_cv_tuneid_"
      outfile <- file.path(opt$res_path, paste0("tuneid_", sprintf("%03d", tune_id), "_fold_", sprintf("%02d", fold_id), ".rds"))
      saveRDS(out, file = outfile)
    }
  } 
  list(fit = fit, sparseFit = sparseFit) 
}

#######################
testSpacemap <- function(test, fit, refit, ...) {
  opt <- list(...)
  resid <- test$Y  - (test$Y%*%fit$Beta + test$X%*%fit$Gamma)
  rssScore <- sum(resid*resid) 
  
  #convergence tolerance already accounted for in refit step
  tol <- if (refit) { opt$tol } else { 0.0 }
  dfParCor <- nonZeroUpper(fit$ParCor, tol)
  dfGamma <- nonZeroWhole(fit$Gamma, tol) 
  degFree <- dfParCor + dfGamma    
  
  matrix(c(rssScore, degFree, dfParCor, dfGamma, fit$deltaMax), 
         nrow = 1, ncol = 5)
}

#######################
testSpace <- function(test, fit, refit, ...) { 
  opt <- list(...)
  if (is.null(opt$residConditional)) opt$residConditional <- TRUE
  #convergence tolerance already accounted for in refit step
  tol <- if (refit) { opt$tol } else { 0.0 }
  if (opt$residConditional)  { 
    resid <- test$Y - (test$Y %*% fit$Beta[opt$Yindex,opt$Yindex] - test$X %*% fit$Beta[opt$Xindex,opt$Yindex])
    dfParCor <- nonZeroUpper(fit$ParCor[opt$Yindex, opt$Yindex], tol)
    dfGamma <- nonZeroWhole(fit$ParCor[opt$Xindex, opt$Yindex], tol) 
    degFree <- nonZeroUpper(fit$ParCor, tol)
  } else {
    resid <- test$XY  - test$XY%*%fit$Beta
    dfParCor <- nonZeroUpper(fit$ParCor, 0.0)
    dfGamma <- 0.0
    degFree <- dfParCor
  }
  rssScore <- sum(resid*resid) 
  
  matrix(c(rssScore, degFree, dfParCor, dfGamma, fit$deltaMax), 
         nrow = 1, ncol = 5)
}

#######################
testModel <- function(test, fit, sparseFit, method, refit, ...) {
  
  opt <- list(...)
  
  if (sparseFit) {
    #currenlty the 0LS refit is returning the partial correlation, not the \Beta's
    fit$Beta <- Beta.coef(fit$ParCor[upper.tri(fit$ParCor)], fit$sig.fit)
    #do not regress on self
    diag(fit$Beta) <- 0 
    if (method == "spacemap") {
      cvScore <- testSpacemap(test = test, fit = fit, refit = refit, ...)
    } else if (method  == "space") {
      cvScore <- testSpace(test = test, fit = fit, refit = refit, ...)
    }
  }
  else {
    cvScore <- matrix(NA, nrow = 1, ncol = 5)
    cvScore[1,5] <- fit$deltaMax
  }
  cvScore
}

dataPart <- function(f, trainIds, testIds, data, 
                      method = c("spacemap", "space")) {
  #define training and test data sets.
  train <- list(); test <- list()
  if (method == "spacemap") {
    train$X <- data$X[trainIds[[f]],]; train$Y <- data$Y[trainIds[[f]],];
    test$X <- data$X[testIds[[f]],]; test$Y <- data$Y[testIds[[f]],];
  } else if (method == "space") { 
    train$XY <- data$XY[trainIds[[f]],]; train$XY <- data$XY[trainIds[[f]],];
    #need this for RSS-score
    test$XY <- data$XY[testIds[[f]],]; test$XY <- data$XY[testIds[[f]],];
    #need this for Conditional Log Likelihood Score
    test$X <- data$X[testIds[[f]],]; test$Y <- data$Y[testIds[[f]],];
  } 
  list(train=train, test=test) 
}

#######################
doCV <- function(data, fold, trainIds, testIds,
                 method = c("spacemap", "space"), tuneGrid, modelFit = FALSE, refit = TRUE, ...) {
  library(foreach)
  
  if (!modelFit) {
    message("Computing CV scores over the grid...")
    cvScores <- foreach(f = seq_len(fold)) %:%
      foreach(l = seq_len(nrow(tuneGrid)), .combine = 'rbind') %dopar% {
        parts <- dataPart(f = f, trainIds = trainIds, testIds = testIds, 
                          data = data, method = method)
        trained <- trainModel(train = parts$train, method = method, tuneGrid = tuneGrid[l,], refit = refit, 
                              fold_id = f, tune_id = l, ...)
        testModel(test = parts$test, fit = trained$fit, sparseFit = trained$sparseFit, method = method, refit = refit, ...) 
      }
  } else {
    message("Computing model fit over CV-selected tuning parameters.")
    cvScores <- foreach(f = seq_len(fold)) %:%
      foreach(l = seq_len(nrow(tuneGrid))) %dopar% {
        parts <- dataPart(f = f, trainIds = trainIds, testIds = testIds, 
                          data = data, method = method)
        trained <- trainModel(train = parts$train, method = method, tuneGrid = tuneGrid[l,], refit = refit,
                              fold_id = f, tune_id = l, ...)
        if(!trained$fit$convergence) { 
          message("CV-selected model did not converge!")
          NA
        } else {
          trained$fit
        }
      }
  }
  cvScores
}

structureScores <- function(cvScores, fold, method) { 

  metrics <- c("rss", "df", "dfParCor", "dfGamma", "deltaMax")
  library(foreach)
  metricScores <- foreach(m = seq_along(metrics)) %do% {
    metricMatrix <- foreach(f = seq_len(fold), .combine = 'cbind') %do% {
      cvScores[[f]][,m]
    } 
  }
  names(metricScores) <- metrics
  metricScores
}

############################
#Find the min score index
#
#@param cvScoresAvg matrix of averaged CV scores. 
#
#@example 
#cvScoresAvg <- matrix(c(1,2,2,1, 11,10,9,9, 3,2,2,5, 8, 7, 9, 10), nrow = 4, ncol = 4)
#colnames(cvScoresAvg) <- c("rss", "bicCondLike", "bicLike", "df")
#minScoreIndex(cvScoresAvg)
minScoreIndex <- function(cvScoresAvg) {
  #========required for more than one score type==========#
  #numScoreType <- ncol(cvScoresAvg) - 1
  #minScores<- apply(cvScoresAvg[,!(colnames(cvScoresAvg) %in% "df")] , 2, min)
  #minIndices <- lapply(1:numScoreType, function(type) which(cvScoresAvg[,type] == minScores[type]))
  #we prefer the more parsimonious model in the event that there are equal scores 
  #minIndex <- sapply(seq_along(minIndices), function(i) minIndices[[i]][which.min(cvScoresAvg[ minIndices[[i]],"df"])])
  #names(minIndex) <- names(minScores)
  #======================================================#
  minRssIds <- which(cvScoresAvg[,"rss"] == min(cvScoresAvg[,"rss"]))
  minRssModelIds <- which.min(cvScoresAvg[minRssIds,"df"])
  c(rss = minRssIds[minRssModelIds])
} 

cvVote <- function(foldFits, method, majorThresh = 0.5, refit, ...) {
  opt <- list(...)
 
  #all below-tolerance and still non-zero coefficients have been zeroed out in trainModel 
  vtol <- 0.0 
  #if (refit) { 
  #  vtol <- 0.0
  #} else {
  #  if(is.null(opt$tol))  {
  #    warning("Did not specify tolerance (tol) as additional argument Defaulting to 1e-6.")
  #    opt$tol <- 1e-6
  #  }
  #  vtol <- opt$tol
  #}
  
  q <- ncol(foldFits[[1]][[1]]$ParCor)
  ParCorVote <- matrix(0, nrow = q, ncol = q)
  if (method ==  "spacemap") {
    p <- nrow(foldFits[[1]][[1]]$Gamma)
    GammaVote <- matrix(0, nrow = p, ncol = q)
  }
   
  for (i in seq_along(foldFits)) {
    if (method == "spacemap") {
      GammaVote <- GammaVote + (abs(foldFits[[i]][[1]]$Gamma) > vtol)
    }
    ParCorVote <- ParCorVote + (abs(foldFits[[i]][[1]]$ParCor) > vtol)
  }
  
  voteFit <- list()
  nconv <- length(foldFits)
  cutoff <- floor(majorThresh*nconv)
  if (method == "spacemap") {
    voteFit$Gamma <- (GammaVote > cutoff) + 0
  }
  voteFit$ParCor <- (ParCorVote > cutoff) + 0
  voteFit$nconv <- nconv
  voteFit$majorThresh <- majorThresh 
  voteFit
}

cv1sd <- function(cvScores, modelSizes) { 
  #1-sd rule
  sdRuleThresh <- mean(cvScores) + sd(cvScores) / sqrt(length(cvScores))
  lessOneSd <- which(cvScores < sdRuleThresh)
  oneSdId <- lessOneSd[which.min(modelSizes[lessOneSd])]
  oneSdId
}


#' Cross validation (CV.Vote) for spacemap and space models.
#'
#' Selects and returns best-tuned model under CV.Vote. 
#'
#' @param data list containing named elements of at least X and Y. Element XY is required when \code{method = "space"}. 
#'X is an (N x P) matrix containing N samples and P predictor variables. 
#'Y is an (N X Q) matrix containing N samples and Q response variables. 
#'XY is an (N X (P + Q)) matrix generated by \code{cbind(X,Y)}.
#' @param fold integer reporting the number of cross validation folds (e.g. 10)
#' @param trainIds list (named optional) of size \code{fold} with the ith element containing 
#'the sample indices of corresponding to the ith training set. 
#' @param testIds list (named optional) of size \code{fold} with the ith element containing 
#'the sample indices of corresponding to the ith test set. 
#' @param method character vector reporting to fit either Spacemap or Space. 
#' @param tuneGrid data.frame named with columns \code{slasso, rlasso, rgroup} when \code{method = "spacemap"}. 
#'Each row in the data.frame corresponds to a tuning parameter set that is input into Spacemap. 
#'So the size of the tuning grid is \code{nrow(tuneGrid)}. When \code{method = "space"}, the previous description still applies, 
#'but is a data.frame with the only column being  \code{slasso}.
#' @param refit logical telling to refit the model after convergence with either ordinary least squares or ridge regression. 
#'Must specify \code{refitRidge = 0.0} as an additional argument to indicate OLS, otherwise\code{refitRidge = 0.001} as a default. 
#' @param ... Additional arguments for Spacemap and Space. For both methods specify values for \code{iter, cd_iter, tol, refitRidge}. 
#'For Space,  \code{sridge, Xindex, Yindex} should be specified.
#'Set \code{sridge = 0.0} unless elastic net penalty is desired.
#'\code{Xindex} are the indices of the X variables in the \code{XY} element of \code{data}. 
#'\code{Yindex} are the indices of the Y variables in the \code{XY} element of \code{data}. 
#'
#' @return List with elements (i) 
#' \code{cvFits} which contains the cv-selected fitted models according to cv-selection type. 
#' \code{tuneGrid} the selected tuning parameter set. 
#' \code{minIndices} the minimum CV-score index of the input tuning grid. 
#' \code{metricScores} A list of Residual sums of squares; total degrees of freedom, Partial Correlation degrees of freedom, 
#' Gamma degrees of freedom, and the final maximum difference between parameter estimates useful for convergence diagnostics.
#'  Each of these is laid out in a matrix of size (tuning grid X number of CV folds). 
#'
#'@examples 
#'
#'\dontrun{
#'
#'
#'#load data
#'data(sim1)
#'dat <- sim1[c("XY", "X", "Y", "Xindex", "Yindex")]
#'N <- nrow(dat$X)
#'P <- ncol(dat$X)
#'Q <- ncol(dat$Y)

#'#set up parallel environment
#'dopar <- TRUE
#'if (dopar) { 
#'  suppressPackageStartupMessages(library(doParallel))
#'  suppressPackageStartupMessages(library(parallel))
#'  ncores <- detectCores()  - 1
#'  cl <- makeCluster(ncores)
#'  registerDoParallel(cl)
#'}
#'
#'#set up training and test sets
#'K <- 10L
#'set.seed(265616L)
#'allSets <- seq_len(N)
#'library(caret)
#'testSets <- caret::createFolds(allSets, k = K)
#'trainSets <- lapply(testSets, function(s) setdiff(allSets, s))
#'
#'#create tuning grid for spacemap
#'tmap <- expand.grid(slasso = c(65,70,73),
#'                    rlasso = c(27,28.8, 30),
#'                    rgroup = c(11, 12.38, 14))
#'#cross validation of spacemap
#'out <- spacemap::crossValidation(dat = dat, fold = K,
#'                                 trainIds = trainSets,  testIds = testSets, 
#'                                 method = "spacemap", tuneGrid = tmap, 
#'                                 iter = 3, tol = 1e-6, cd_iter = 1e7, 
#'                                 sridge = 0, refitRidge = 0)
#'}
crossValidation <- function(data, fold, trainIds, testIds,
                            method = c("spacemap", "space"), tuneGrid, refit = TRUE, ...) {
  
  ## compute cross-validated metric scores for each fold and for each tune set
  foldScores <- doCV(data = data, fold = fold, trainIds = trainIds, 
                   testIds = testIds, method = method, 
                   tuneGrid = tuneGrid, modelFit = FALSE, refit = refit, ...)
  
  ##restructure list$fold[tune index, metric index ] to be list$metric$[tune index, fold index ] 
  metricScores <- structureScores(cvScores = foldScores, fold = fold, method = method)
  
  ##Average cross-validated scores across folds
  #while accounting for non-convergent folds
  testSetLen <- sapply(testIds, length)
  total <- sum(testSetLen)
  foldconvmat <- !is.na(metricScores$rss)
  nfoldconv <- rowSums(foldconvmat)
  ntestconv <- apply(X = foldconvmat, MARGIN  = 1, FUN = function(not_na) sum(testSetLen[not_na]))
  metricScoresAvg <- foreach(i = seq_along(metricScores), .combine = 'cbind') %do% {
    score <- metricScores[[i]]
    stype <- names(metricScores)[[i]]
    msum <- apply(X = score , MARGIN = 1, FUN = sum, na.rm = TRUE) 
    
    ###Selection of tuning grid restricted to case when majority of test hold outs
    ###were evaluated (chunked by folds) 
    
    # divide rss by number of test sets with corresponding trained convergence
    #divide other metrics by the number of converged folds
    
    if (stype == "rss") { 
      ifelse(ntestconv > 0.5*total, msum / ntestconv, Inf)
    } else {
      ifelse(ntestconv > 0.5*total, msum / nfoldconv, Inf)
    }
  }
  colnames(metricScoresAvg) <- names(metricScores)
    
  ##Find the minimizing tuning parameter set
  minIndices <- minScoreIndex(metricScoresAvg)
  if (metricScoresAvg[minIndices[1],"rss"] == Inf) { 
    stop("No convergence for any fold or any tuning parameter selection.")
  }   

  tuneGridParamLabels <- names(tuneGrid)
  #get fitted models according to best score for a given score-type.
  cvFits <- foreach(j = seq_along(minIndices)) %do% {
    #make sure subset of data.frame is still a data.frame (for case where there is one column)
    bestTune <- as.data.frame(tuneGrid[minIndices[j],])
    names(bestTune) <- tuneGridParamLabels

    #only retrain on sets that converged. 
    rss_cv_scores <- metricScores[["rss"]][minIndices[j],]
    converged <- !is.na(rss_cv_scores) 
    nconv <- sum(converged)
    df_cv_scores <- metricScores[["df"]][minIndices[j],]

    ## Return best-fit models: Note add parameter for returning model instead of scores
    foldFits <- doCV(data = data, fold = nconv, trainIds = trainIds[converged], 
                         testIds = testIds[converged], method = method, 
                         tuneGrid = bestTune, modelFit = TRUE, refit = refit, ...)
    if (any(is.na(foldFits))) { 
      warning("Retraining models: initial convergence did not converge second time! Do not trust cvVote, cv1sd output.")
    }
    #cv 1-sd
    cv1sdId <- cv1sd(cvScores = rss_cv_scores[converged], 
                     modelSizes = df_cv_scores[converged])
    cv1sdFit <- foldFits[[cv1sdId]][[1]]
    #cv vote
    voteFit <- cvVote(foldFits = foldFits, method = method, refit = refit, ...)
    #K fold cross validation (now train on full data).
    cvK <- trainModel(train = data, method = method, tuneGrid = tuneGrid[minIndices[j],], 
                      refit = refit, fold_id = 0 , tune_id = 0, ...)$fit
    list(cvK = cvK, vote = voteFit, cv1sd = cv1sdFit)
  }
  names(cvFits) <- names(minIndices)
  list(cvFits = cvFits, tuneGrid = tuneGrid[minIndices,], minIndices = minIndices, metricScores = metricScores)
}

cvPerf <- function(cvOut, trueParCor, method, conditional = TRUE, tol = 1e-6, Xindex=NULL, Yindex=NULL) {
  library(foreach)
  cvTypes <- names(cvOut$cvFits$rss)
  perfOut <- foreach(metric = names(cvOut$cvFits), .combine = 'rbind') %do% {
    metricStep <- foreach(cvType = names(cvOut$cvFits$rss), .combine = 'rbind') %do% {
      if (method == "spacemap") {
        perf <- spacemap::reportJointPerf(fitParCor = list(xy = cvOut$cvFits[[metric]][[cvType]]$Gamma, 
                                                           yy = cvOut$cvFits[[metric]][[cvType]]$ParCor),
                                          trueParCor = trueParCor, tol = tol, 
                                          verbose = FALSE)  
      } else if (method == "space") {
        if (conditional) { 
          perf <- spacemap::reportJointPerf(fitParCor = list(xy = cvOut$cvFits[[metric]][[cvType]]$ParCor[Xindex,Yindex], 
                                                             yy = cvOut$cvFits[[metric]][[cvType]]$ParCor[Yindex,Yindex]),
                                            trueParCor = trueParCor, tol = tol, 
                                            verbose = FALSE)
        } else { 
          perf <- spacemap::reportPerf(cvOut$cvFits[[metric]][[cvType]]$ParCor,
                                       trueParCor$yy, YY = TRUE, tol = tol, 
                                       verbose = FALSE)
        }
      }
      perf
    }
    metricStep <- as.data.frame(metricStep)
    cbind(metricStep, metric = metric, cvType = cvTypes)
  }
  rownames(perfOut) <- NULL
  if (conditional) { 
    perfOut[,c(1:7, 20, 21)]
  } else { 
    perfOut
  }
}
