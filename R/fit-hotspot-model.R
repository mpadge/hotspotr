#' fit_hotspot_model
#'
#' Fits a neutral model to an observed series of hotspot values in terms of
#' standard deviation of environmental variables, spatial autocorrelation
#' parameters, and number of iterations of spatial autocorrelation.
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
#' @return A vector of three values as estimated by the neutral model:
#' \enumerate{
#'   \item sd0 = standard deviation of normal distribution
#'   \item alpha = temporal autocorrelation coefficient
#'   \item niters = number of iterations of spatial autocorrelation
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

    maxnbs <- max (sapply (nbs, length))
    nbsi <- lapply (seq (maxnbs), function (i) get_nbsi (i, nbs))

    # Get reference distributions of raw values and associated AC coefficients
    zs <- sort (z, decreasing=TRUE)
    zs <- (zs - min (zs)) / diff (range (zs))
    acs <- rcpp_ac_stats (z, nbs, wts, ac_type)
    acs <- (acs - min (acs)) / diff (range (acs))

    opt_fn <- function (x)
    {
        # x [1] = alpha
        # x [2] = sd0
        # x [3] = log_scale
        # x [4] = niters
        # x [5] = ntests
        dat <- rcpp_neutral_hotspots_ntests (nbs, wts, nbsi, 
                                             alpha=x [1], sd0=x [2],
                                             log_scale=x[3], niters=x [4], 
                                             ac_type=ac_type, ntests=x [5])
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
    log_scale <- TRUE # constant throughout
    search_span <- 20 # length of vectors used for search
    niter_lims <- c (1, 10)
    alpha_lims <- c (-1, 1)
    sd_lims <- c (1e-6, 1)
    niters <- unique (round (seq (niter_lims [1], niter_lims [2], 
                                  length.out=search_span)))
    alpha <- seq (alpha_lims [1], alpha_lims [2], length.out=search_span)
    sd0 <- seq (sd_lims [1], sd_lims [2], length.out=search_span)
    tol <- 1e-3
    trial <- 1
    maxtrials <- 10
    niters1 <- round (mean (niters))
    alpha1 <- mean (alpha) # arbitrary values used to calculate err
    sd1 <- mean (sd0)
    err <- ymin_old <- opt_fn (c (alpha1, sd1, log_scale, niters, ntests))
    while (err > tol & trial <= maxtrials)
    {
        if (verbose) message ('trial#', trial, ': ', appendLF=FALSE)
        niters2 <- niters1
        alpha2 <- alpha1
        sd2 <- sd1

        yniters <- sapply (niters, function (i) 
                           opt_fn (c (alpha1, sd1, log_scale, i, ntests)))
        if (all (diff (yniters) > 0) | all (diff (yniters) < 0) |
            length (niters) < 6)
        {
            niters1 <- niters [which.min (yniters)]
        } else
        {
            mod <- loess (yniters ~ niters, span=0.75, degree=1)$fitted
            niters1 <- niters [which.min (mod)]
        }

        yac <- sapply (alpha, function (i) 
                       opt_fn (c (i, sd1, log_scale, niters1, ntests)))
        mod <- loess (yac ~ alpha, span=0.5)$fitted
        alpha1 <- alpha [which.min (mod)]

        ysd <- sapply (sd0, function (i) 
                       opt_fn (c (alpha1, i, log_scale, niters1, ntests)))
        mod <- loess (ysd ~ sd0, span=0.5)$fitted
        sd1 <- sd0 [which.min (mod)]

        niter_lims [1] <- round ((niter_lims [1] + niters1) / 2)
        niter_lims [2] <- round ((niter_lims [2] + niters1) / 2)
        alpha_lims [1] <- (alpha_lims [1] + alpha1) / 2
        alpha_lims [2] <- (alpha_lims [2] + alpha1) / 2
        sd_lims [1] <- (sd_lims [1] + sd1) / 2
        sd_lims [2] <- (sd_lims [2] + sd1) / 2
        niters <- unique (round (seq (niter_lims [1], niter_lims [2], 
                                      length.out=search_span)))
        alpha <- seq (alpha_lims [1], alpha_lims [2], length.out=50)
        sd0 <- seq (sd_lims [1], sd_lims [2], length.out=50)

        ymin <- min (mod)
        err <- abs (alpha1 - alpha2) + abs (sd1 - sd2)
        if (ymin > ymin_old)
        {
            message ('stopping search because error is increasing')
            trial <- maxtrials
        }
        err <- min (c (err, ymin_old - ymin))
        ymin_old <- ymin
        trial <- trial + 1
        if (verbose)
            message ('(niters, alpha, sd) = (', niters1, ', ', alpha1, ', ', 
                     sd1, '), err = ', err)
    }

    #list (sd0=op$par [1], ac=op$par [2])
    list (sd0=sd1, ac=alpha1, niters=niters1)
}
