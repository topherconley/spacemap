#pass cross-validation defaults if not specified by user
variTrainParam <- function(...) {
  opt <- list(...)
  if (is.null(opt$iter)) { opt$iter <- 3L } 
  if (is.null(opt$tol)) { opt$tol <- 1e-6 } 
  if (is.null(opt$cdmax)) { opt$cdmax <- 1e7 } 
  if (is.null(opt$sridge)) { opt$sridge <- 0.0 }
  if (is.null(opt$weight)) { opt$weight <- NULL }
  if (is.null(opt$remWeight)) { opt$remWeight <- NULL }
  if (is.null(opt$sig)) { opt$sig <- NULL }
  if (is.null(opt$rho)) { opt$rho <- NULL }
  #hedge against multicolinearity, for numerical stability
  if (is.null(opt$refitRidge)) {  opt$refitRidge <- 0.01 } 
  opt
}

#######################
trainModel <- function(train, method, tuneGrid, refit = TRUE, fold_id, tune_id, resPath, 
                       givenX, Xindex, Yindex, aszero, ...) {
  
  
  opt <- variTrainParam(...)
  
  #feedback on ols refit assumption of model sparsity
  sparseFit <- TRUE
  requireNamespace("Matrix")
  if (method == "spacemap") {
    fit <- spacemap(Y = train$Y, X = train$X, 
                    lam1 = tuneGrid$lam1, sridge = opt$sridge, 
                    lam2 = tuneGrid$lam2, lam3 = tuneGrid$lam3, 
                    sig=opt$sig, weight=opt$weight, remWeight = opt$remWeight, iter= opt$iter,
                    tol = opt$tol, cdmax = opt$cdmax, iscale = FALSE)
    #zero out those below tolerance. 
    fit$ParCor[abs(fit$ParCor) <= aszero ] <- 0.0
    fit$Gamma[abs(fit$Gamma) <= aszero ] <- 0.0
    #reject models that did not converge
    if(!fit$convergence) return(list(fit = fit, sparseFit = FALSE))
    #recompute ols estimates based on training set: note change fit object by reference
    if(refit) {
      sparseFit <- olsRefitSpacemap(Y = train$Y, X = train$X, 
                                    ParCor = fit$ParCor, Gamma = fit$Gamma, 
                                    sigma = fit$sig.fit, RSS = fit$rss, tol = opt$tol, 
                                    ridge = opt$refitRidge)
    }

    #write to file
    out <- list(ThetaXY = Matrix(fit$Gamma), ThetaYY = Matrix(fit$ParCor), sig = fit$sig.fit,
                tune_id = tune_id, fold_id = fold_id, convergence = fit$convergence, sparseFit = sparseFit)
    outfile <- file.path(resPath, paste0("tuneid_", sprintf("%03d", tune_id), 
                                          "_fold_", sprintf("%02d", fold_id), 
                                          ".rds"))
    saveRDS(out, file = outfile)
    
  } else if (method == "space") {
    fit <- spacemap::space(Y = train$XY, lam1 = tuneGrid, sridge = opt$sridge, iter = opt$iter, 
                       cdmax = opt$cdmax, tol = opt$tol, iscale = FALSE,
                       sig = opt$sig, rho = opt$rho)
    #zero out those below tolerance. 
    fit$ParCor[abs(fit$ParCor) <= aszero] <- 0.0
    if(!fit$convergence) return(list(fit = fit, sparseFit = FALSE))
    if(refit) {
      sparseFit <- olsRefitSpace(Y = train$XY,
                                 ParCor = fit$ParCor, 
                                 sigma = fit$sig.fit, RSS = fit$rss, tol = opt$tol,
                                 ridge = opt$refitRidge)
    }
  
    #write to file
    out <- list(sig = fit$sig.fit,
                tune_id = tune_id, fold_id = fold_id, 
                convergence = fit$convergence, sparseFit = sparseFit)
    
    if (givenX) { 
      xy <- fit$ParCor[Xindex,Yindex]
      yy <- fit$ParCor[Yindex,Yindex]
      out$ThetaXY = Matrix(xy) 
      out$ThetaYY = Matrix(yy)
    } else {
      out$ThetaYY = Matrix(fit$ParCor) 
    }
    outfile <- file.path(resPath, paste0("tuneid_", sprintf("%03d", tune_id), 
                                          "_fold_", sprintf("%02d", fold_id), 
                                          ".rds"))
    saveRDS(out, file = outfile)
  } 
  list(fit = fit, sparseFit = sparseFit) 
}

#######################
testSpacemap <- function(test, fit, refit) {
  resid <- test$Y  - (test$Y%*%fit$Beta + test$X%*%fit$Gamma)
  rssScore <- sum(resid*resid) 
  #convergence tolerance already accounted for train step
  tol <- 0.0
  dfParCor <- nonZeroUpper(fit$ParCor, tol)
  dfGamma <- nonZeroWhole(fit$Gamma, tol) 
  degFree <- dfParCor + dfGamma    
  matrix(c(rssScore, degFree, dfParCor, dfGamma, fit$deltaMax), 
         nrow = 1, ncol = 5)
}

#######################
testSpace <- function(test, fit, refit, givenX, Xindex, Yindex) { 
  #convergence tolerance already accounted for in train step
  tol <- 0.0
  if (givenX)  { 
    resid <- test$Y - (test$Y %*% fit$Beta[Yindex,Yindex] - test$X %*% fit$Beta[Xindex,Yindex])
    dfParCor <- nonZeroUpper(fit$ParCor[Yindex, Yindex], tol)
    dfGamma <- nonZeroWhole(fit$ParCor[Xindex, Yindex], tol) 
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
testModel <- function(test, fit, sparseFit, method, refit, givenX, Xindex, Yindex) {
  
  if (sparseFit) {
    #currenlty the 0LS refit is returning the partial correlation, not the \Beta's
    fit$Beta <- Beta.coef(fit$ParCor[upper.tri(fit$ParCor)], fit$sig.fit)
    #do not regress on self
    diag(fit$Beta) <- 0 
    if (method == "spacemap") {
      cvScore <- testSpacemap(test = test, fit = fit, refit = refit)
    } else if (method  == "space") {
      cvScore <- testSpace(test = test, fit = fit, refit = refit,
                           givenX = givenX, Xindex = Xindex, Yindex = Yindex)
    }
  }
  else {
    cvScore <- matrix(NA, nrow = 1, ncol = 5)
    cvScore[1,5] <- fit$deltaMax
  }
  cvScore
}

dataPart <- function(f, trainIds, testIds, data, 
                      method = c("spacemap", "space"), givenX) {
  #define training and test data sets.
  train <- list(); test <- list();
  if (method == "spacemap") {
    train$X <- data$X[trainIds[[f]],]; train$Y <- data$Y[trainIds[[f]],];
    test$X <- data$X[testIds[[f]],]; test$Y <- data$Y[testIds[[f]],];
  } else if (method == "space") { 
    train$XY <- data$XY[trainIds[[f]],];
    if(givenX) { 
      test$X <- data$X[testIds[[f]],]; test$Y <- data$Y[testIds[[f]],];  
    } else { 
      test$XY <- data$XY[testIds[[f]],];
    }
  }
  list(train=train, test=test) 
}


structureScores <- function(cvScores, fold, method) { 

  metrics <- c("rss", "df", "dfParCor", "dfGamma", "deltaMax")
  requireNamespace("foreach")
  #for R CMD check NOTE passing
  m <- NULL; f <- NULL;
  metricScores <- foreach(m = seq_along(metrics)) %do% {
    metricMatrix <- foreach(f = seq_len(fold), .combine = 'cbind') %do% {
      cvScores[[f]][,m]
    } 
  }
  names(metricScores) <- metrics
  metricScores
}

averageScores <- function(metricScores, testIds) { 
  ##Average cross-validated scores across folds
  #while accounting for non-convergent folds
  testSetLen <- sapply(testIds, length)
  total <- sum(testSetLen)
  foldconvmat <- !is.na(metricScores$rss)
  nfoldconv <- rowSums(foldconvmat)
  ntestconv <- apply(X = foldconvmat, MARGIN  = 1, FUN = function(not_na) sum(testSetLen[not_na]))
  i <- 0L #for R CMD check 
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
  metricScoresAvg
}

############################
# Find the min score index
minScoreIndex <- function(cvScoresAvg) {
  #if there are multiple equally good scores, take the most sparse model.
  minRssIds <- which(cvScoresAvg[,"rss"] == min(cvScoresAvg[,"rss"]))
  minRssModelIds <- which.min(cvScoresAvg[minRssIds,"df"])
  c(rss = minRssIds[minRssModelIds])
} 


#' Cross validation (CV.Vote) for spacemap and space models.
#'
#' Selects and returns best-tuned model under CV.Vote. 
#' 
#' @inheritParams spacemap
#' @param trainIds List of integer vectors, where each integer vector contains 
#' a split of training sample indices pertaining to \code{Y, X}. 
#' @param testIds List of integer vectors, where each integer vector contains 
#' a split of  test sample indices pertaining to \code{Y, X}. Required to 
#' be of the same length as \code{trainIds}. 
#' @param method Character vector indicates network inference with function 
#' \code{\link{spacemap}} when \code{method = "spacemap"} or function 
#' \code{\link{space}} when \code{method = "space"}. If \code{X} is 
#' non-null and \code{method = "space"}, then \code{space} will 
#' infer (x--x, x--y, y--y) edges but only report (x--y, y--y) edges. 
#' @param tuneGrid Named with columns \code{lam1, lam2, lam3} 
#' when \code{method = "spacemap"}. Each row in the data.frame corresponds to a
#'  tuning parameter set that is input into \code{\link{spacemap}}. 
#' When \code{method = "space"}, supply a data.frame with only one column 
#' being \code{lam1}.
#' @param resPath Character vector specifying the directory where each 
#' model fit is written to file through serialization by \code{saveRDS}. 
#' Defaults to temporary directory that will be deleted at end of the R session. 
#' It is recommended to specify a directory where results can be stored permanently. 
#' @param refit Logical indicates to refit the model after convergence to
#' reduce bias induced by penalty terms. Default to TRUE. The refit step
#' defaults to a  ridge regression with small penalty of 0.01 to 
#' encourage numerical stability. The user can change the ridge
#' penalty by adding an additional parameter \code{refitRidge}. 
#' @param thresh Numerical threshold between 0 and 1 (defaults to 0.5 or majority vote) 
#' indicating the minimum proportion of times () a given edge must be represented 
#' in the trained models to be reported in the final CV.Vote model. For example, 
#' If 0.5 is specified, and there are 10 training splits, then an edge in the final 
#' model must be reported in 6 of the 10 traiing models. 
#' @param aszero Positive numeric value (defaults to 1e-6) indicating at what point to consider extremely 
#' small parameter estimates of \eqn{\Gamma} and \eqn{\rho} as zero. 
#' @param ... Additional arguments for \code{\link{spacemap}} or \code{\link{space}}
#' to change their default settings (e.g. setting \code{tol = 1e-4}). 
#'
#' @return  A list containing  
#'  \itemize{  
#'  \item A list called \code{cvVote} with two elements: 
#'  \itemize{
#'   \item \code{xy}, an adjacency matrix where \eqn{xy(p,q)} element
#'   is 1 for an edge between \eqn{x_p} and \eqn{y_q} and 0 otherwise; and
#'   \item \code{yy} Adjacency matrix where \eqn{yy(q,l)} element
#'   is 1 for an edge between \eqn{y_q} and \eqn{y_l} and 0 otherwise. 
#'   }
#'  \item \code{minTune} List containing the optimal tuning penalty set. 
#'  \item \code{minIndex} Integer specifying the index of \code{minTune} in \code{tuneGrid}. 
#'  \item \code{metricScores} Data.frame for input to \code{\link{tuneVis}} for
#'  inspecting the CV score curve and model size as a function of the tuning penalties.
#' }
#' @import RcppArmadillo
#' @import Rcpp
#' @export
#' @examples 
#'library(spacemap)
#'data(sim1)
#'##########################
#'#DEFINE TRAINING/TEST SETS
#'library(caret)
#'#sample size
#'N <- nrow(sim1$X)
#'#number of folds
#'K <- 5L
#'set.seed(265616L)
#'#no special population structure, but create randomized dummy structure of A and B
#'testSets <- createFolds(y = sample(x = c("A", "B"), size = N, replace = TRUE), k = K)
#'trainSets <- lapply(testSets, function(s) setdiff(seq_len(N), s))
#'nsplits <- sapply(testSets, length)
#'##########################
#'#SPACE (Y input)
#'tsp <- expand.grid(lam1 = seq(65, 75, length = 3))
#'cvspace <- cvVote(Y = sim1$Y, 
#'                  trainIds = trainSets, testIds = testSets, 
#'                  method = "space", tuneGrid = tsp) 
#'##########################
#'# SPACEMAP (Y and X input)
#'tmap <- expand.grid(lam1 = seq(65, 75, length = 2), 
#'                    lam2 = seq(21, 35, length = 2), 
#'                    lam3 = seq(10, 40, length = 2))
#'cvsmap <- cvVote(Y = sim1$Y, X = sim1$X, 
#'                 trainIds = trainSets, testIds = testSets, 
#'                 method = "spacemap", tuneGrid = tmap)
#'                 
cvVote <- function(Y, X = NULL, trainIds, testIds,
                            method = c("spacemap", "space"), tuneGrid,  
                            resPath = tempdir(), refit = TRUE, 
                            thresh = 0.5,  iscale = TRUE, aszero = 1e-6, ...) {
  ################
  #check arguments
  ################
  #paired test/train splits
  stopifnot(length(trainIds) == length(testIds))
  fold <- length(trainIds)
  method <- match.arg(method)
  
  if(!is.matrix(Y)) stop("Y is not a matrix.")
  givenX <- !is.null(X) 
  Xindex <- NULL 
  Yindex <- NULL
  if(method == "space" & givenX)  {
    if(!is.matrix(X)) stop("X is not a matrix.")
    if(iscale) { 
      X <- scale(X)
      Y <- scale(Y)
    }
    data <- list(XY = cbind(X, Y), Y = Y, X = X)
    Xindex <- seq_len(ncol(X))
    Yindex <- (ncol(X) + 1):(ncol(X) + ncol(Y))
  } else if(method == "space" & !givenX){ 
    if(iscale) {
      Y <- scale(Y)
    }
    data <- list(XY = Y)
  } else if (method == "spacemap") { 
    if(!is.matrix(X)) stop("X is not a matrix.")
    if(iscale) { 
      X <- scale(X)
      Y <- scale(Y)
    }
    data <- list(Y = Y, X = X)
  }
  
  ## compute cross-validated metric scores for each fold and for each tune set
  requireNamespace("foreach")
  message("Computing CV scores over the grid...")
  iscale <- FALSE
  #for R CMD check NOTE passing
  f <- NULL; l <- NULL;
  foldScores <- foreach(f = seq_len(fold)) %:% 
    foreach(l = seq_len(nrow(tuneGrid)), .combine = 'rbind', .packages = "spacemap") %dopar% {
      parts <- dataPart(f = f, trainIds = trainIds, testIds = testIds, 
                        data = data, method = method, givenX = givenX)
      trained <- trainModel(train = parts$train, method = method, tuneGrid = tuneGrid[l,], refit = refit, 
                            fold_id = f, tune_id = l, resPath = resPath,
                            givenX = givenX, Xindex = Xindex, Yindex = Yindex, aszero = aszero, ...)
      tmp <- testModel(test = parts$test, fit = trained$fit, sparseFit = trained$sparseFit, 
                       method = method, refit = refit, 
                       givenX = givenX, Xindex = Xindex, Yindex = Yindex) 
      tmp
    }
  
  ##restructure list$fold[tune index, metric index ] to be list$metric$[tune index, fold index ] 
  metricScores <- structureScores(cvScores = foldScores, fold = fold, method = method)
  
  #average across folds
  metricScoresAvg <- averageScores(metricScores = metricScores, testIds = testIds)
  
  ##Find the minimizing tuning parameter set
  minIndex <- minScoreIndex(metricScoresAvg)
  if (metricScoresAvg[minIndex[1],"rss"] == Inf) { 
    stop("No convergence for any fold or any tuning parameter selection.")
  }   
  
  #obtain model fit files from best tune index
  files <- list.files(path = resPath, 
                           pattern =  paste0("tuneid_", sprintf("%03d", minIndex)),
                           full.names = TRUE)
  voteFit <- cvVoteRDS(files = files, tol = 0.0, 
                       thresh = thresh, method = method,
                       givenX = givenX)
  
  if (method == "spacemap") { 
    mtg <- tuneGrid[minIndex,]
    minTune  <- list(lam1 = mtg$lam1, lam2 = mtg$lam2, lam3 = mtg$lam3)  
  } else if (method == "space") { 
    minTune  <- list(lam1 = tuneGrid[minIndex,])  
  }
  
  list(cvVote = voteFit, minTune = minTune, minIndex = minIndex, 
       logcvScore = log(metricScoresAvg[minIndex,"rss"]),
       metricScores = metricScores, bestFoldFiles = files)
}

##################
#PERFORMANCE
##################

#' CV.Vote Model Performance
#' 
#' Convenience function for evaluating edge detection performance of 
#' CV.Vote model fits. 
#' 
#' @inheritParams reportJointPerf
#' @param cvOut Output from \code{\link{cvVote}}.
#' @param method Character vector indicates \code{cvOut} came from
#' spaceMap or SPACE model. See parameter \code{method} in \code{\link{cvVote}}.
#' @param givenX Logical (defaults to FALSE). Set to TRUE for \code{method == "space"}
#' when \code{cvOut} had non-null \code{X} input to \code{\link{cvVote}}.
#' @param Xindex Integer vector of indices of X variables in partial correlation matrix. Defaults to 
#' NULL, but must be non-null when \code{method == "space" & givenX = TRUE}. 
#' @param Yindex Integer vector of indices of Y variables in partial correlation  matrix. Defaults to 
#' NULL, but must be non-null when \code{method == "space" & givenX = TRUE}. 
#' @seealso \code{\link{reportJointPerf}}, \code{\link{reportPerf}}
#' @export
cvPerf <- function(cvOut, truth, method, givenX = FALSE, Xindex=NULL, Yindex=NULL) {
  if (method == "spacemap") {
    perf <- spacemap::reportJointPerf(fit = list(xy = cvOut$cvVote$xy, 
                                                       yy = cvOut$cvVote$yy),
                                      truth = truth, tol = 1e-6, 
                                      verbose = FALSE)  
  } else if (method == "space") {
    if (givenX) { 
      if(is.null(Xindex) | is.null(Yindex)) { 
        stop("Must specify Xindex and Yindex; cannot be NULL.")  
      }
      perf <- spacemap::reportJointPerf(fit = list(xy = cvOut$cvVote$yy[Xindex,Yindex], 
                                                         yy = cvOut$cvVote$yy[Yindex,Yindex]),
                                        truth = truth, tol = 1e-6, 
                                        verbose = FALSE)
    } else { 
      perf <- spacemap::reportPerf(cvOut$cvVote$yy,
                                   truth$yy, YY = TRUE, tol = 1e-6, 
                                   verbose = FALSE)
    }
  }
  perf
}
