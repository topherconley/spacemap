

setwd("~/repos/bioc_spacemap/")
library(devtools)

devtools::check()
devtools::clean_dll()
devtools::install()


library(lintr)
lintr::lint_package()





