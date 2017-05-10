
#' Visualize cross validation metrics
#' 
#' Diagnose the quality of the tuning grid used in \code{\link{cvVote}}
#' through plotting  CV scores and network size against the
#' tuning penalties. 
#' 
#' @param cvOut Output from \code{\link{cvVote}}. 
#' @param testSetLen Integer vector specifying length of each test set associated with 
#' \code{cvOut}. 
#' @param tuneParam1 Numeric vector specifying the tuning penalty placed on the 
#' x-axis. 
#' @param tuneParam2 Numeric vector specifying the tuning penalty placed on the 
#' x-axis. Defaults to NULL. If non-null, the output is a heatmap where 
#' CV scores are colored intensities on a 2-D grid of \code{tuneParam1,tuneParam2}. 
#' @param tuneParam1Name Character vector labelling the X-axis with the specified tuning parameter. 
#' @return If \code{tuneParam2 == NULL}
#' Return a list of 4 \code{ggplot2} objects. 
#' \enumerate{
#' \item Scatter plot of log (CV score) vs. \code{tuneParam1}. 
#' \item Scatter plot of "Avg. Total # of Edges" vs. \code{tuneParam1}. 
#' \item Scatter plot of "Avg. # of Edges (Y -- Y)" vs. \code{tuneParam1}.
#' \item Scatter plot of "Avg. # of Edges (x -> Y)" vs. \code{tuneParam1}.
#' }
#' 
#' If \code{tuneParam2 != NULL}
#' Returns a heatmap (see \code{ggplot2::geom_tile}) where 
#' CV scores are colored intensities on a 2-D grid of \code{tuneParam1,tuneParam2}.
#' @importFrom ggplot2 qplot ggplot aes geom_tile theme_bw xlab ylab
#' @importFrom reshape2 melt
#' @importFrom foreach %do%
#' @export
#' @seealso \code{\link{cvVote}}
tuneVis <- function(cvOut, testSetLen,
                    tuneParam1, tuneParam1Name = c("lam1", "lam2", "lam3"),
                    tuneParam2 = NULL) {
  
  requireNamespace("ggplot2")
  tuneParam1Name <- match.arg(tuneParam1Name)
  total <- sum(testSetLen)
  foldconvmat <- !is.na(cvOut$metricScores[["rss"]])
  nfoldconv <- rowSums(foldconvmat)
  ntestconv <- apply(X = foldconvmat, MARGIN  = 1, FUN = function(not_na) sum(testSetLen[not_na]))
  requireNamespace("foreach")
  #for R CMD check NOTE passing
  i <- NULL;
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
  
  if (is.null(tuneParam2)) {
    dat <- as.data.frame(avgScores)
    g1 <- qplot(data = dat, x = tuneParam1, y = log(rss), geom = "point") + 
      ylab("log (CV Score)") + xlab(tuneParam1Name) + theme_bw()
    g2 <- qplot(data = dat, x = tuneParam1, y = df, geom = "point") +
      ylab("Avg. Total # of Edges") + xlab(tuneParam1Name) + theme_bw()
    g3 <- qplot(data = dat, x = tuneParam1, y = dfParCor, geom = "point") + 
      ylab("Avg. # of Edges (Y -- Y)") + xlab(tuneParam1Name) + theme_bw()
    g4 <- qplot(data = dat, x = tuneParam1, y = dfGamma, geom = "point") + 
      ylab("Avg. # of Edges (X -> Y)") + xlab(tuneParam1Name) + theme_bw()
    #requireNamespace("gridExtra")
    #gg1 <- grid.arrange(g1,g2,g3,g4, ncol = 2)
    gg1 <- list(g1,g2,g3,g4)
#     diagnose <- as.data.frame(cbind(tuneName = tuneParam1, avgScores))
#     requireNamespace("reshape2")
#     ggdiagnose <- melt(diagnose, id.vars = "tuneName", 
#                        measure.vars = colnames(avgScores), value.name = "avg_score",
#                        variable.name  = "score_type")
#     
#     requireNamespace("ggplot2")
#     gg1 <- ggplot(data = ggdiagnose, aes(x = tuneName, y = avg_score, colour = score_type)) + geom_point()  + facet_grid(. ~ score_type)
#     #gg2 <- qplot(tuneParam, rowMeans(cvOut$metricScores$df))
  } else {
    diagnose <- as.data.frame(cbind(tuneName1 = tuneParam1, tuneName2 = tuneParam2, rss = avgScores[,"rss"]))
    requireNamespace("reshape2")
    ggdiagnose <- melt(diagnose, id.vars = c("tuneName1", "tuneName2"), 
                       measure.vars = "rss", value.name = "rss")
    gg1 <- ggplot(data = ggdiagnose, aes(x = tuneName1, y = tuneName2, fill = rss)) + geom_tile()
  }
  #list(gg1, gg2)
  gg1
}
