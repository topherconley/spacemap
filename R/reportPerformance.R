#'@title Report Model Performance
#'@description Assesses the performance of space, remMap, and spacemap models 
#'against a known truth by reporting power and false discovery rate. 
#'@param fitParCor A numeric matrix encoding the fitted model parameter estimates. 
#'@param trueParCor A numeric matrix encoding the true model parameters.
#'@param  YY A logical value indicating that \code{fitParCor} estimates
#' have the same structure as the \emph{Space} model. If \code{YY=TRUE}, then 
#' the fitted and the truth inputs must be symmetric matrices. 
#'@param tol A numeric value specifying a lower bound on the tolerance of a non-zero 
#'parameter. 
#'@param verbose A logical value indicating whether to report the performance to the console. 
#'@return A numeric, named vector where: (X,Y) is the joint edge performance of X-Y and Y-Y edges. 
#'\itemize{
#'\item{power}{Joint power on (X,Y)}
#'\item{fdr}{Joint false discovery rate (FDR) on (X,Y)}
#'\item{powerXY}{Power of X-Y edges}
#'\item{fdrXY}{FDR of X-Y edges}
#'\item{powerYY}{Power of Y-Y edges}
#'\item{fdrYY}{FDR of Y-Y edges}
#'\item{tp}{True positives of (X,Y)}
#'\item{fn}{False negatives of (X,Y)}
#'\item{fp}{False positives of (X,Y)}
#'\item{tpXY}{True positives of X-Y}
#'\item{fnXY}{False negatives of X-Y}
#'\item{fpXY}{False positives of X-Y}
#'\item{tpYY}{True positives of Y-Y}
#'\item{fnYY}{False negatives of Y-Y}
#'\item{fpYY}{False positives of Y-Y}
#'}
#'
reportPerf <- function(fitParCor, trueParCor, YY = TRUE, tol = 1e-6, 
                       verbose = TRUE) {
  
  #true adjacency
  trueAdj <- abs(trueParCor) > tol
  #fitted adjacency 
  fitAdj <- abs(fitParCor) > tol
  #verify that the dimensions between comparisons match
  stopifnot(nrow(fitAdj) == nrow(trueAdj), ncol(fitAdj) == ncol(trueAdj))
  
  if (YY) {
    #do not count self-adjacency
    diag(trueAdj) <- FALSE
    diag(fitAdj) <- FALSE
    
    #required to be square and symmetric.
    stopifnot(nrow(fitAdj) == ncol(fitAdj), 
              nrow(trueAdj) == ncol(trueAdj),
              all(trueAdj == t(trueAdj)),
              all(fitAdj == t(fitAdj)))
    
    #only care about the upper triangle (symmetric)
    ufit <- fitAdj[upper.tri(fitAdj)]
    utrue <- trueAdj[upper.tri(trueAdj)]
    
    #true positives
    tp <- as.integer(sum(ufit == T & utrue == T))
    #true negatives
    tn <- as.integer(sum(ufit == F & utrue == F))
    #false positives
    fp <- as.integer(sum(ufit == T & utrue == F))
    #false negatives
    fn <- as.integer(sum(ufit == F & utrue == T))
  } else {
    tp <- as.integer(sum(fitAdj == T & trueAdj == T))
    tn <- as.integer(sum(fitAdj == F & trueAdj == F))
    fp <- as.integer(sum(fitAdj == T & trueAdj == F))
    fn <- as.integer(sum(fitAdj == F & trueAdj == T))
  }

  #sensitivity = power
  power <- algoPower(tp, fn)
  #sensitivity (not used, overoptimistic in sparse settings because tn so big)
  #spec <- tn / (tn + fp)
  fdr <- algoFDR(tp, fp)
  #mathews correlation coefficient, not sensitive to prevalence
  mcc <- algoMCC(tp, fn, fp, tn)
  
  
  if (verbose)
    message("Power: ", round(power,4), "\tFDR: ", round(fdr,4))
  invisible((c(power=power, fdr=fdr, mcc=mcc, tp=tp, fn=fn, fp=fp, tn=tn)))
}

reportJointPerf <- function(fitParCor, trueParCor, tol = 1e-6, 
                       verbose = TRUE) {
  #report performance metrics for XY, YY results respectively.
  xy <- reportPerf(fitParCor$xy, trueParCor$xy, YY = FALSE, tol = tol, 
                               verbose = FALSE)
  yy <- reportPerf(fitParCor$yy, trueParCor$yy, YY = TRUE, tol = tol, 
                   verbose = FALSE)
  
  #combine xy and yy results into one power and fdr calculation
  tp <-  xy["tp"] + yy["tp"] 
  fn <-  xy["fn"] + yy["fn"] 
  fp <-  xy["fp"] + yy["fp"]
  tn <-  xy["tn"] + yy["tn"]
  
  power <- algoPower(tp, fn)
  fdr <- algoFDR(tp, fp)
  mcc <- algoMCC(tp, fn, fp, tn)
  
  if (verbose)
    message("Power: ", round(power,4), "\tFDR: ", round(fdr,4))
  perf <- c(power, fdr, mcc,
            xy["power"], xy["fdr"], 
            yy["power"], yy["fdr"],
            yy["mcc"], xy["mcc"],
            tp, fn, fp, tn,
            xy["tp"], xy["fn"], xy["fp"],
            yy["tp"], yy["fn"], yy["fp"])
  names(perf) <- c("power", "fdr", "mcc",
                   "powerXY", "fdrXY", 
                   "powerYY", "fdrYY", 
                   "mccYY", "mccXY",
                   "tp", "fn", "fp", "tn",
                   "tpXY", "fnXY", "fpXY",
                   "tpYY", "fnYY", "fpYY")
  invisible(perf)
}

#'@title Statistical Power of an Algorithm
#'@param tp The number of true positive classifications.
#'@param fn The number of false negative classifications.
algoPower <- function(tp , fn) {
  tp / (tp + fn)
}

#'@title False Discovery Rate of an Algorithm
#'@param tp The number of true positive classifications.
#'@param fp The number of false positive classificationsm.
algoFDR <- function(tp, fp) {
  if ((tp + fp) == 0) 
    fdr <- 0
  else 
    fdr <- fp / (tp + fp)
  return(fdr)
}

#'@title Mathew's Correlation Coefficient of Classifier
#'@param tp The number of true positive classifications.
#'@param fn The number of false negative classifications.
#'@param fp The number of false positive classifications.
#'@param tn The number of true negative classifications.
algoMCC <- function(tp, fn, fp, tn) {
  #make sure is integer
  #stopifnot(is.integer(tp),
  #          is.integer(fn),
  #          is.integer(fp),
  #          is.integer(tn))
  #make sure is non-negative
  stopifnot(c(tp,fn,tp,tn) > -1)

  numer <- tp*tn - fp*fn
  denomDefault <- c(tp + fp, tp + fn, tn + fp, tn + fn)
  if (all(denomDefault > 0)) { 
    mcc <- numer / sqrt(prod(denomDefault))
  } else {
    mcc <- numer
  }
  return(mcc)
}
