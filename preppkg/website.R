
setwd("~/repos/spacemap/")
#install.packages(c("devtools", "ggplot2", "Rcpp", 
#                  "RcppArmadillo", "foreach", "doRNG", 
#                   "Matrix", "igraph", "crayon"))
library(devtools)
load_all()
document()
install()
devtools::clean_dll()
#install_github(repo = "hadley/pkgdown", force = T)
#install.packages("rmarkdown")
library(pkgdown)
library(spacemap)
build_site(preview = T)
#build_home()
#build_reference()


build(binary = T, vignettes = F)
?build


install.packages(c("Rcpp", "RcppArmadillo", 
                   "foreach", "doRNG", 
                   "Matrix", "igraph"))


### build on windows
"C:\Program Files\R\R-3.3.2\bin\R.exe CMD INSTALL --build compile-both spacemap/"
remove.packages("spacemap")
install.packages(pkgs = "C:/Users/topher/Shared/repos/spacemap_0.45.0.zip", 
                 repos = NULL)

##build on linux
setwd("~/repos/spacemap/")
library(devtools)
build(binary = T, vignettes = F)


## build on mac

library(spacemap)

