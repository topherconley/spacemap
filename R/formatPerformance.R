formatPerf <- function(fit, truth, method = "spacemap", aszero = 1e-6) { 
  comparison <- factor(c("(yy,xy)", "yy", "xy"), levels = c( "(yy,xy)", "yy", "xy"))
  #performance names to index
  pn <- list(power = c("power", "powerYY", "powerXY"),
             fdr = c("fdr", "fdrYY", "fdrXY"), 
             mcc = c("mcc", "mccYY", "mccXY"),
             tp = c("tp", "tpYY", "tpXY"),
             fn = c("fn", "fnYY", "fnXY"),
             fp = c("fp", "fpYY", "fpXY"))
  library(spacemap)
  perf <- reportJointPerf(fit, truth, aszero = aszero, verbose = FALSE)
  data.frame(method = method,
             power = perf[pn$power], fdr =  perf[pn$fdr], mcc = perf[pn$mcc],
             tp = perf[pn$tp], fn = perf[pn$fn], fp = perf[pn$fp],
             comparison = comparison)
}
