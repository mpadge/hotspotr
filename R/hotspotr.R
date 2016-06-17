#' hotspotr
#'
#' Analyses whether the statistical properties of a spatial pattern of hotspots
#' may be reproduced with a simple neutral model. 
#'
#' @section Functions:
#' \tabular{ll}{
#' \code{gearyc}\tab Geary's C statistic on square grid\cr
#' \code{getisord}\tab Getis-Ord spatial association statistic on square grid\cr
#' \code{ives2d}\tab Simulate 2D square grid using model of Ives & Klopfer
#' (Ecology 1997)\cr
#' \code{morani}\tab Moran's I statistic on square grid\cr
#' \code{neutral1d}\tab Neutral model in one dimension\cr
#' \code{neutral2d}\tab Neutral model in two dimensions\cr
#' \code{run_tests}\tab Test observed data with range of possible (1-D & 2-D)
#' models\cr
#' \code{test1d}\tab Test observed data against a one dimensional neutral
#' model\cr
#' \code{test2d}\tab Test observed data against a two dimensional neutral
#' model\cr
#' }
#'
#' @name hotspotr
#' @docType package
#' @importFrom stats rnorm t.test wilcox.test optim optimise
#' @importFrom graphics lines legend title par plot.new
#' @importFrom msm rtnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib hotspotr
NULL
