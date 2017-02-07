avgMetrics <- function(cvOut, testSetLen) {
  total <- sum(testSetLen)
  foldconvmat <- !is.na(cvOut$metricScores[["rss"]])
  nfoldconv <- rowSums(foldconvmat)
  ntestconv <- apply(X = foldconvmat, MARGIN  = 1, FUN = function(not_na) sum(testSetLen[not_na]))
  library(foreach)
  avgScores <- foreach(i = seq_along(cvOut$metricScores), .combine = 'cbind') %do% {
    score <- cvOut$metricScores[[i]]
    stype <- names(cvOut$metricScores)[[i]]
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
  colnames(avgScores) <- names(cvOut$metricScores)
  avgScores
}
