#for devtools and pkdown functions
setwd("~/repos/spacemap/")

#############################
#UPDATE R PACKAGE
#############################
library(devtools)
#load_all()
#document()
#install(quick = TRUE, dependencies = T)
#devtools::clean_dll()

#############################
#BUILD WEBSITE 
#############################
#install_github(repo = "hadley/pkgdown", force = T)
library(pkgdown)
library(spacemap)
build_site(preview = T)
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

### build on windows
#build_win()

#run this in command prompt
"C:\Program Files\R\R-3.3.3\bin\R.exe CMD INSTALL --build compile-both spacemap/"
remove.packages("spacemap")
install.packages(pkgs = "C:/Users/topher/Shared/repos/spacemap_0.51.0.zip", 
                 repos = NULL)

##build on Linux an MAC OS
setwd("~/repos/spacemap/")
library(devtools)
build(binary = T, vignettes = F)
