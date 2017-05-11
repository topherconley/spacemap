#' Network toolkit example input data 
#'
#' 
#' 
#' @description Input for illustrating how to apply spaceMap's 
#' network toolkit to data from the 
#' Breast Cancer Proteogenomics Landscape Study (BCPLS)
#' by Mertins et. al., (2016, Nature). 
#'  
#' @details The list \code{bcpls} is described in detail at
#'  \url{https://topherconley.github.io/spacemap/articles/neta.html#input}.
#'
#' 
#' @docType data
#' @keywords datasets
#' @name bcpls
#' @usage data(bcpls)
#' @format A list containing objects:
#' \itemize{
#'  \item \code{net} spaceMap's Boot.Vote ensemble network for the BCPLS data.
#'  \item \code{yinfo} Annotations for proteins (response variables).
#'  \item \code{xinfo} Annotations for genomic copy number alterations (predictor variables).
#'  \item \code{bdeg} Degree distributions of bootstrap ensemble replicates.
#'  \item \code{go2eg} Mapping from Gene Ontology (GO) to proteins. 
#' }
NULL