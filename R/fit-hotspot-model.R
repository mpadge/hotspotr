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
#' xy <- cbind (rep (seq (size), each=size), rep (seq (size), size))
#' dhi <- 1 # for rook; dhi=1.5 for queen
#' nbs <- spdep::dnearneigh (xy, 0, dhi)
#' z <- runif (length (nbs))
#' test <- fit_hotspot_model (z=z, nbs=dat$nbs, alpha=0.1, sd=0.1, ntests=100)
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

    # Get reference distributions of raw values and associated AC coefficients
    zs <- sort (z, decreasing=TRUE)
    zs <- (zs - min (zs)) / diff (range (zs))
    acs <- rcpp_ac_stats (z, nbs, wts, ac_type)
    acs <- (acs - min (acs)) / diff (range (acs))

    opt_fn <- function (x)
    {
        # x [1] = alpha
        # x [2] = sd0
        # x [3] = niters
        # x [4] = ntests
        dat <- neutral_hotspots_ntests2 (nbs=nbs, wts=wts, alpha=x[1], sd0=x[2],
                                         niters=x[3], ntests=x[4])
        sum ((dat [,1] - zs) ^ 2) + sum ((dat [,2] - acs) ^ 2)
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
    tol <- 1e-3
    iter <- 0
    maxiters <- 10
    alpha1 <- mean (alpha) # arbitrary values used to calculate err
    sd1 <- mean (sd0)
    niters <- 1
    err <- opt_fn (c (alpha1, sd1, niters, ntests))
    while (err > tol & iter < maxiters)
    {
        if (verbose) message ('iteration#', iter, ': ', appendLF=FALSE)
        alpha2 <- alpha1
        sd2 <- sd1

        yac <- sapply (alpha, function (i) opt_fn (c (i, sd1, niters, ntests)))
        mod <- loess (yac ~ alpha, span=0.5)$fitted
        alpha1 <- alpha [which.min (mod)]
        ysd <- sapply (sd0, function (i) opt_fn (c (alpha1, i, niters, ntests)))
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

    #list (sd0=op$par [1], ac=op$par [2])
    list (sd0=sd1, ac=alpha1)
}
