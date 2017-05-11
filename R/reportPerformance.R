#'@title Report performance of x-y or y-y edge detection
#'@description Assesses the performance of space and spaceMap
#'against a known truth by reporting power and false discovery rate
#'for either x-y OR y-y edge detection. 
#'@param fit A numeric matrix encoding the fitted model parameter estimates. 
#'@param truth A numeric matrix encoding the true model parameters.
#'@param  YY A logical defaulting to TRUE indicating that 
#'\code{fit} and \code{truth} are  symmetric matrices 
#' corresponding to estimated and true partial correlations for y-y edges. 
#' Otherwise, assume \code{fit} and \code{truth} corresponds to the estimated and true
#' \eqn{\Gamma} regression coefficient matrix for x-y edges.
#'@param aszero  A numeric value specifying the point at which a parameter estimate should 
#'be effectively considered of zero value. 
#'@param verbose A logical value indicating whether to report the performance to the console. 
#'@export
#'@return A numeric, named vector containing
#'\itemize{
#' \item power or sensitivity
#' \item fdr or false discovery rate
#' \item mcc or Matthew's Correlation Coefficient
#' \item tp  or true positive
#' \item fn or false negative
#' \item fp or false positive
#' \item tn or true negative
#'}
#'@seealso reportJointPerf
reportPerf <- function(fit, truth, YY = TRUE, aszero = 1e-6, 
                       verbose = TRUE) {
  
  #true adjacency
  trueAdj <- abs(truth) > aszero
  #fitted adjacency 
  fitAdj <- abs(fit) > aszero
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

#'@title Report joint performance of (x-y,y-y) edge detection
#'@description Assesses the network learning performance of spaceMap (or space)
#'against a known truth by reporting power and false discovery rate. 
#'@param fit A list of numeric matrices named 'xy' and 'yy' 
#' encoding the fitted model parameter estimates (or the true adjacency matrices). 
#'@param truth A list of numeric matrices named 'xy' and 'yy' 
#' encoding the true model parameters (or an adjacency matrix).
#'@param aszero A numeric value specifying the point at which a parameter estimate should 
#'be effectively considered of zero value. 
#'@param verbose A logical value indicating whether to report the performance to the console. 
#'@seealso reportPerf
#'@export
#'@return A numeric, named vector 
#'\itemize{
#'\item{power}{ Joint power on (x-y,y-y)}
#'\item{fdr}{ Joint false discovery rate (FDR) on (x-y,y-y)}
#'\item{powerXY}{ Power of x-y edges}
#'\item{fdrXY}{ FDR of x-y edges}
#'\item{powerYY}{ Power of y-y edges}
#'\item{fdrYY}{ FDR of y-y edges}
#'\item{tp}{ True positives of (x-y,y-y)}
#'\item{fn}{ False negatives of (x-y,y-y)}
#'\item{fp}{ False positives of (x-y,y-y)}
#'\item{tpXY}{ True positives of x-y}
#'\item{fnXY}{ False negatives of x-y}
#'\item{fpXY}{ False positives of x-y}
#'\item{tpYY}{ True positives of y-y}
#'\item{fnYY}{ False negatives of y-y}
#'\item{fpYY}{ False positives of y-y}
#'}
#'
reportJointPerf <- function(fit, truth, aszero = 1e-6, 
                       verbose = TRUE) {
  #report performance metrics for XY, YY results respectively.
  xy <- reportPerf(fit$xy, truth$xy, YY = FALSE, aszero = aszero, 
                               verbose = FALSE)
  yy <- reportPerf(fit$yy, truth$yy, YY = TRUE, aszero = aszero, 
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


#' Power function
#' 
#' Compute the power given true positives and false negatives.
#' 
#'@param tp The number of true positive classifications.
#'@param fn The number of false negative classifications.
algoPower <- function(tp , fn) {
  tp / (tp + fn)
}

#' False discovery rate
#' 
#' Compute the FDR given true positives and false positives.
#' 
#'@param tp The number of true positive classifications.
#'@param fp The number of false positive classificationsm.
algoFDR <- function(tp, fp) {
  if ((tp + fp) == 0) 
    fdr <- 0
  else 
    fdr <- fp / (tp + fp)
  return(fdr)
}

#' Matthew's Correlation Coefficient (MCC)
#' 
#' Compute the MCC given the confusion matrix.
#' 
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
