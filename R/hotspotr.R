#' hotspotr
#'
#' Analyses whether the statistical properties of a spatial pattern of hotspots
#' may be reproduced with a simple neutral model. 
#'
#' @section Functions:
#' \tabular{ll}{
#' \code{ives2d}\tab Simulate model of Ives & Klopfer (Ecology 1997)\cr
#' \code{neutral2d}\tab Neutral model in two dimensions\cr
#' \code{run_tests}\tab Test observed data with range of possible (1-D & 2-D)
#' models\cr
#' \code{test2d}\tab Test observed data against a two dimensional neutral
#' model\cr
#' }
#'
#' @name hotspotr
#' @docType package
#' @importFrom stats rnorm t.test wilcox.test optim optimise
#' @importFrom methods is
#' @importFrom graphics lines legend title par plot.new text
#' @importFrom msm rtnorm
#' @importFrom spdep dnearneigh
#' @importFrom Rcpp evalCpp
#' @useDynLib hotspotr
NULL
