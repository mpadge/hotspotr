#' fit_hotspot_model
#'
#' Fits a neutral model to an observed series of hotspot values in terms of
#' temporal and spatial autocorrelation parameters, and number of iterations of
#' these.
#'
#' @param z Vector of observed values to be tested
#' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
#' point
#' @param wts Weighting factors for each neighbour; must have same length as
#' nbs. Uniform weights used if not given.
#' @param alpha Vector of two components providing starting values for the
#' strength of autocorrelation in time and space
#' @param ntests Number of repeats of neutral model used to calculate mean
#' rank--scale distribution
#' @param ac_type type of autocorrelation statistic to use in tests
#' (\code{moran}, \code{geary}, or \code{getis-ord}=\code{go})
#' @param sd0 Standard deviation of truncated normal distribution used to model
#' environmental variation (with mean of 1)
#' @param verbose If TRUE, dump progress details to screen
#' @param plot If TRUE, produces a plot of rank--scale distributions
#'
#' @return A vector of four values as estimated by the neutral model:
#' \enumerate{
#'   \item alpha respectively containing temporal and spatial autocorrelation
#'   coefficients
#'   \item nt = Number of successive layers of spatio-temporal autocorrelation
#'   required to reproduce statistical properties of observed data
#' }
#'
#' @section Note:
#' Fitting these neutral models is **not** a standard optimisation problem
#' because the models are very noisy. Although \code{optim} with
#' \code{method="SANN"} may be used, it often generates extremely large values
#' for \code{alpha} (for example, > 10). \code{DEoptim} could also be applied,
#' yet in generally does not explore anything useful---if given starting
#' parameters, it will generally remain exactly in that place.
#'
#' The approach employed here reflects the comment of
#' https://stat.ethz.ch/pipermail/r-help/2015-May/428751.html
#' through simply producing regular series, fitting loess models, and taking the
#' corresponding minima.
#'
#' @examples
#' \dontrun{
#' alpha <- c (0.1, 0.1)
#' dat <- ives (size=10, nt=10, sd0=0.1, alpha=alpha)
#' test <- fit_hotspot_model (z=dat$dat$z, nbs=dat$nbs, alpha=alpha, ntests=100)
#' }
#'
#' @export
fit_hotspot_model <- function (z, nbs, wts, alpha=c(0.1, 0.1), ntests=100,
                               ac_type='moran', sd0=0.1, verbose=FALSE, plot=FALSE)
{
    if (!is.numeric (z)) 
        stop ('z must be numeric')
    if (!is (nbs, 'nb'))
        stop ('nbs must of class spdep::nb')
    if (length (z) != length (nbs))
        stop ('nbs must have same length as z')

    ac_type <- tolower (ac_type)
    if (substring (ac_type, 1, 1) == 'g')
    {
        if (substring (ac_type, 3, 3) == 'a')
            ac_type <- 'geary'
        else
            ac_type <- 'getis_ord'
    } else
        ac_type <- 'moran'

    size <- length (z)

    if (missing (wts)) 
        wts <- lapply (nbs, function (x) rep (1, length (x)) / length (x))

    ac <- rcpp_ac_stats (z, nbs, wts, ac_type)
    zs <- sort ((z - min (z)) / diff (range (z)), decreasing=TRUE)
    test <- NULL # remove no visible binding warning

    # Initial 3D optimisation to get nt. This is a repeat of the code from
    # neutral-hotspots-ntests, repeated here to avoid constantly recreating the
    # cluster for each call
    clust <- parallel::makeCluster (parallel::detectCores () - 1)
    parallel::clusterExport (clust, list ("nbs", "wts", "alpha", "sd0", 
                                          "ac_type"), envir=environment ())
    invisible (parallel::clusterCall (clust, function () {
                            while (length (grep ('hotspotr', getwd ())) > 0) 
                                setwd ("..")
                            devtools::load_all ("hotspotr")
                            setwd ("./hotspotr")
                                             }))

    fn_n <- function (x)
    {
        parallel::clusterExport (clust, list ("x"), envir=environment ())
        z <- parallel::parLapply (clust, seq (ntests), function (i) 
                        {
                            z1 <- rcpp_neutral_hotspots (nbs, wts, 
                                                         alpha_t=alpha [1],
                                                         alpha_s=alpha [2],
                                                         sd0=sd0, nt=x)
                            ac1 <- rcpp_ac_stats (z1, nbs, wts, ac_type)
                            z1 <- (sort (z1, decreasing=TRUE) - min (z1)) / 
                                    diff (range (z1)) 
                            rbind (z1, ac1)
                        })
        ac1 <- colMeans (do.call (rbind, lapply (z, function (i) i [2,])))
        z1 <- colMeans (do.call (rbind, lapply (z, function (i) i [1,])))
        sum ((z1 - zs) ^ 2) + sum ((ac1 - ac) ^ 2)
    }

    if (verbose) cat ("optimising for n ... ")
    nt <- 1:20
    nt0 <- 1
    while (nt0 == 1 | nt0 == 20)
    {
        y <- sapply (nt, fn_n)
        mod <- loess (y ~ nt, span=0.5)$fitted
        nt0 <- nt [which.min (mod)]
    }
    nt <- nt0
    parallel::stopCluster (clust)

    # then reduce to 2d optimisation
    clust <- parallel::makeCluster (parallel::detectCores () - 1)
    parallel::clusterExport (clust, list ("nbs", "wts", "sd0", "nt",
                                          "ac_type"), envir=environment ())
    invisible (parallel::clusterCall (clust, function () {
                            while (length (grep ('hotspotr', getwd ())) > 0) 
                                setwd ("..")
                            devtools::load_all ("hotspotr")
                            setwd ("./hotspotr")
                                             }))
    fn_a <- function (x)
    {
        parallel::clusterExport (clust, list ("x"), envir=environment ())
        z <- parallel::parLapply (clust, seq (ntests), function (i) 
                        {
                            z1 <- rcpp_neutral_hotspots (nbs, wts, 
                                                         alpha_t=x [1],
                                                         alpha_s=x [2],
                                                         sd0=sd0, nt=nt)
                            ac1 <- rcpp_ac_stats (z1, nbs, wts, ac_type)
                            z1 <- (sort (z1, decreasing=TRUE) - min (z1)) / 
                                    diff (range (z1)) 
                            rbind (z1, ac1)
                        })
        ac1 <- colMeans (do.call (rbind, lapply (z, function (i) i [2,])))
        z1 <- colMeans (do.call (rbind, lapply (z, function (i) i [1,])))
        sum ((z1 - zs) ^ 2) + sum ((ac1 - ac) ^ 2)
    }
    if (verbose) cat ("done.\noptimising for alpha ... ")
    # Note that a full 2D search with loess surface is much more likely to to
    # find minima at the edges rather than the interior. While nevertheless
    # notionally better, it's much quicker to conduct a series of 1D searches
    at <- -10:10 / 10
    at <- at [!at == 0]
    asp <- at
    alpha_s <- 0.1
    for (i in 1:4)
    {
        y <- sapply (at, function (i) fn_a (c (i, alpha_s)))
        mod <- loess (y ~ at, span=0.5)$fitted
        alpha_t <- at [which.min (mod)]

        y <- sapply (asp, function (i) fn_a (c (alpha_t, i)))
        mod <- loess (y ~ asp, span=0.5)$fitted
        alpha_s <- asp [which.min (mod)]

        if (i == 2)
        {
            at <- alpha_t + -10:10 / 50
            at <- at [!at == 0]
            asp <- alpha_s + -10:10 / 50
            asp <- asp [!asp == 0]
        }
    }

    if (verbose) cat ("done\n")

    parallel::stopCluster (clust)

    list (alpha=c(alpha_t, alpha_s), nt=nt)
}
