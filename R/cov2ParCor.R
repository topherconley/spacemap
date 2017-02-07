

cov2ParCor <- function(Sigma) {
  #validate  argument is a matrix and symmetry.
  stopifnot(is.matrix(Sigma), identical(Sigma, t(Sigma)))
  SigmaInv <- solve(Sigma)
  #see src/ functions of Rcpp. 
  conc2parcor(SigmaInv)
}