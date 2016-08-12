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
#' @param ac_type type of autocorrelation statistic to use in tests
#' (\code{moran}, \code{geary}, or \code{getis-ord}=\code{go})
#' @param ntests Number of repeats of neutral model used to calculate mean
#' rank--scale distribution
#' @param verbose If TRUE, dump progress details to screen
#' @param plot If TRUE, produces a plot of rank--scale distributions
#'
#' @return A vector of two values as estimated by the neutral model:
#' \enumerate{
#'   \item sd0 = standard deviation of normal distribution
#'   \item alpha = temporal autocorrelation coefficient
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
fit_hotspot_model <- function (z, nbs, wts, ac_type='moran', ntests=100,
                               verbose=FALSE, plot=FALSE)
{
    if (missing (z)) stop ('z must be provided')
    if (missing (nbs)) stop ('nbs must be provided')
    if (!is.numeric (z)) stop ('z must be numeric')
    if (!is (nbs, 'nb')) stop ('nbs must of class spdep::nb')
    if (length (z) != length (nbs))
        stop ('z must have the same length as nbs')
    if (!missing (wts))
    {
        if (length (wts) != length (nbs))
            stop ('wts must have the same length as nbs')
    } else
        wts <- lapply (nbs, function (x) rep (1, length (x)) / length (x))

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

    get_nbsi <- function (i)
    {
        res <- lapply (seq (nbs), function (j)
                       {
                           if (length (nbs [[j]]) >= i)
                               c (j, nbs [[j]] [i], length (nbs [[j]]))
                           else
                               NULL
                       })
        res <- res [lapply (res, length) != 0]
        res <- do.call (rbind, res)
        data.frame (to=res [,1], from=res [,2], n=res [,3])
    }
    maxnbs <- max (sapply (nbs, length))

    ac <- rcpp_ac_stats (z, nbs, wts, ac_type)
    zs <- sort ((z - min (z)) / diff (range (z)), decreasing=TRUE)

    # Optimise over 2D space of (sd0, alpha)
    clust <- parallel::makeCluster (parallel::detectCores () - 1)
    exports <- list ("size", "nbs", "wts", "ac_type", "alpha", "sd0",
                     "ntests", "get_nbsi", "maxnbs")
    parallel::clusterExport (clust, exports, envir=environment ())
    invisible (parallel::clusterCall (clust, function () {
                            while (length (grep ('hotspotr', getwd ())) > 0) 
                                setwd ("..")
                            devtools::load_all ("hotspotr")
                            setwd ("./hotspotr")
                                             }))

    opt_fn <- function (x)
    {
        parallel::clusterExport (clust, list ("x"), envir=environment ())
        z <- parallel::parLapply (clust, seq (ntests), function (i) 
                                  {
                                      z1 <- msm::rtnorm (size, mean=1, sd=x [1], 
                                                         lower=0, upper=2)

                                      for (j in seq (x [3]))
                                      {
                                          z2 <- rep (0, size ^ 2)
                                          for (k in seq (maxnbs))
                                          {
                                              nbsi <- get_nbsi (k)
                                              z2 [nbsi$to] <- z2 [nbsi$to] + 
                                                  ((1 - x [2]) * z1 [nbsi$to] +
                                                   x [2] * z1 [nbsi$from]) / 
                                                  nbsi$n
                                          }
                                          z1 <- z2
                                      }

                                      ac1 <- rcpp_ac_stats (z1, nbs, wts, ac_type)
                                      z1 <- (sort (z1, decreasing=TRUE) - 
                                             min (z1)) / diff (range (z1)) 
                                      #z1 <- sort (log10 (z1), decreasing=TRUE)
                                      #z1 <- (z1 - min (z1)) / diff (range (z1))
                                      rbind (z1, ac1)
                                  })
        ac1 <- colMeans (do.call (rbind, lapply (z, function (i) i [2,])))
        z1 <- colMeans (do.call (rbind, lapply (z, function (i) i [1,])))
        sum ((z1 - zs) ^ 2) + sum ((ac1 - ac) ^ 2)
    }

    # opt_fn is too noisy for optim to work
    #control <- list (reltol=1e-3)
    #if (verbose)
    #{
    #    control <- c (control, trace=1)
    #    message ('optimising ...')
    #}
    #op <- optim (c (0.1, 0.1), opt_fn, control=control)

    # Instead, a succession of localised 1D loess-smoothed approaches are used
    alpha_lims <- c (-1, 1)
    sd_lims <- c (1e-6, 1)
    alpha <- seq (alpha_lims [1], alpha_lims [2], length.out=50)
    sd0 <- seq (sd_lims [1], sd_lims [2], length.out=50)
    err <- opt_fn (c (mean (sd0), mean (alpha)))
    tol <- 1e-3
    iter <- 0
    maxiters <- 10
    alpha1 <- mean (alpha) # arbitrary values used to calculate err
    sd1 <- mean (sd0)
    while (err > tol & iter < maxiters)
    {
        if (verbose) message ('iteration#', iter, ': ', appendLF=FALSE)
        alpha2 <- alpha1
        sd2 <- sd1

        yac <- sapply (alpha, function (i) opt_fn (c (sd1, i)))
        mod <- loess (yac ~ alpha, span=0.5)$fitted
        alpha1 <- alpha [which.min (mod)]
        ysd <- sapply (sd0, function (i) opt_fn (c (i, alpha1)))
        mod <- loess (ysd ~ sd0, span=0.5)$fitted
        sd1 <- sd0 [which.min (mod)]

        alpha_lims [1] <- (alpha_lims [1] + alpha1) / 2
        alpha_lims [2] <- (alpha_lims [2] + alpha1) / 2
        sd_lims [1] <- (sd_lims [1] + sd1) / 2
        sd_lims [2] <- (sd_lims [2] + sd1) / 2
        alpha <- seq (alpha_lims [1], alpha_lims [2], length.out=50)
        sd0 <- seq (sd_lims [1], sd_lims [2], length.out=50)

        err <- abs (alpha1 - alpha2) + abs (sd1 - sd2)
        iter <- iter + 1
        if (verbose)
            message ('(alpha, sd) = (', alpha1, ', ', sd1, '), err = ', err)
    }

    parallel::stopCluster (clust)

    #list (sd0=op$par [1], ac=op$par [2])
    list (sd0=sd1, ac=alpha1)
}
