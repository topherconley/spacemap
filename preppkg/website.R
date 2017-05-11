#for devtools and pkdown functions
setwd("~/repos/spacemap/")

#############################
#UPDATE R PACKAGE
#############################
library(devtools)
load_all()
document()

library(roxygen2)
roxygenize(clean = TRUE)
install(quick = TRUE, dependencies = T)
devtools::clean_dll()
devtools::check()
#############################
#BUILD WEBSITE 
#############################
#install_github(repo = "hadley/pkgdown", force = T)
library(pkgdown)
library(spacemap)
ptm <- proc.time()
build_site(preview = T)
proc.time() - ptm
build_home()
build_reference()
build_articles()
build_news()

#############################
#BUILD PKG BINARY RELEASE
#############################
build(binary = T, vignettes = F)

#dependencies of spacemap
install.packages(c("Rcpp", "RcppArmadillo", 
                   "foreach", "doRNG", 
                   "Matrix", "igraph"))

########
#WINDOWS
########

#STEP 1: remove dynamically linked libraries
devtools::clean_dll()
#STEP 2: remove any installed version of spacemap,
#and make sure dependencies are installed.
remove.packages("spacemap")


cranpkgs <- c("Rcpp","RcppArmadillo",
             "foreach","Matrix",
             "rngtools", "ggplot2",
             "reshape2")
install.packages(pkgs = cranpkgs)

## try http:// if https:// URLs are not supported
source("http://bioconductor.org/biocLite.R")
biocpkgs <- c("GenomicRanges")
biocLite(biocpkgs)

#STEP 3: build binary
# Run verbatim as administrator in Windows Command Prompt
# "C:\Program Files\R\R-3.3.3\bin\R.exe" CMD INSTALL --build --compile-both spacemap/
#STEP 4: check that binary installs.
install.packages(pkgs = "C:/Users/topher/Shared/repos/spacemap_0.52.0.zip", 
                 lib = "C:/Users/topher/Shared/R/win-library/3.3",
                 repos = NULL)
library(spacemap)

##build on Linux an MAC OS
setwd("~/repos/spacemap/")
library(devtools)
build(binary = T, vignettes = F)
