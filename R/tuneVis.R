

tuneVis <- function(cvOut, testSetLen,
                    tuneParam1, tuneParam2 = NULL, tuneParam1Name = "slasso") {
  
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
  
  if (is.null(tuneParam2)) {
    library(ggplot2)
    dat <- as.data.frame(avgScores)
    g1 <- qplot(data = dat, x = tuneParam1, y = log(rss), geom = "point") + 
      ylab("log (RSS CV Score )") + xlab(tuneParam1Name) + theme_bw()
    g2 <- qplot(data = dat, x = tuneParam1, y = df, geom = "point") +
      ylab("Total Deg. of Freedom") + xlab(tuneParam1Name) + theme_bw()
    g3 <- qplot(data = dat, x = tuneParam1, y = dfParCor, geom = "point") + 
      ylab("Deg. of Freedom (Y -- Y)") + xlab(tuneParam1Name) + theme_bw()
    g4 <- qplot(data = dat, x = tuneParam1, y = dfGamma, geom = "point") + 
      ylab("Deg. of Freedom (X -> Y)") + xlab(tuneParam1Name) + theme_bw()
    library(gridExtra)
    #gg1 <- grid.arrange(g1,g2,g3,g4, ncol = 2)
    gg1 <- list(g1,g2,g3,g4)
#     diagnose <- as.data.frame(cbind(tuneName = tuneParam1, avgScores))
#     library(reshape2)
#     ggdiagnose <- melt(diagnose, id.vars = "tuneName", 
#                        measure.vars = colnames(avgScores), value.name = "avg_score",
#                        variable.name  = "score_type")
#     
#     library(ggplot2)
#     gg1 <- ggplot(data = ggdiagnose, aes(x = tuneName, y = avg_score, colour = score_type)) + geom_point()  + facet_grid(. ~ score_type)
#     #gg2 <- qplot(tuneParam, rowMeans(cvOut$metricScores$df))
  } else {
    diagnose <- as.data.frame(cbind(tuneName1 = tuneParam1, tuneName2 = tuneParam2, rss = avgScores[,"rss"]))
    library(reshape2)
    ggdiagnose <- melt(diagnose, id.vars = c("tuneName1", "tuneName2"), 
                       measure.vars = "rss", value.name = "rss")
    gg1 <- ggplot(data = ggdiagnose, aes(x = tuneName1, y = tuneName2, fill = rss)) + geom_tile()
  }
  #list(gg1, gg2)
  gg1
}
