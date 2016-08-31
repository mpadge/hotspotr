#' neutral_hotspots
#'
#' Implements neutral model of hotspot values and associated autocorrelation
#' statustics. Current version is a simply internal R loop, while
#' \code{neutral_hotspots_ntests} is a parallel version of exactly the same
#' thing.
#'
#' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
#' point
#' @param wts Weighting factors for each neighbour; must have same length as
#' nbs. Uniform weights used if not given.
#' @param alpha strength of spatial autocorrelation 
#' @param sd0 Standard deviation of truncated normal distribution used to model
#' environmental variation (with mean of 1)
#' @param ac_type Type of autocorrelation statistic to use in tests
#' (\code{moran}, \code{geary}, or \code{getis-org}=\code{go})
#' @param niters Number of successive layers of spatial autocorrelation
#' @param log_scale If TRUE, raw hotspot values are log-transformed
#' @param ntests Number of tests over which to generate an average result
#' @param parallel If true, the tests are conducted using the \code{R} package
#' \code{parallel}, otherwise using (non-parallel) code{Rcpp} loops.
#' @param seed Random seed
#'
#' @return A vector of hotspot values sorted from high to low
#'
#' @seealso \code{ives}
#'
#' @examples
#' # First set up a grid of rectangular neighbours
#' size <- 10
#' xy <- cbind (rep (seq (size), each=size), rep (seq (size), size))
#' dhi <- 1 # for rook; dhi=1.5 for queen
#' nbs <- spdep::dnearneigh (xy, 0, dhi)
#' dat <- neutral_hotspots (nbs, ntests=1000)
#'
#' @export
neutral_hotspots <- function (nbs, wts, alpha=0.1, sd0=0.1, ac_type='moran', 
                              niters=1, log_scale=TRUE, ntests=100,
                              parallel=FALSE, seed)
{
    if (missing (nbs)) stop ('nbs must be given')

    if (missing (wts)) 
        wts <- lapply (nbs, function (x) rep (1, length (x)) / length (x))

    #if (alpha [1] <= 0)
    #    stop ('neutral model only makes sense with finite temporal autocorrelation')

    ac_type <- tolower (ac_type)
    if (substring (ac_type, 1, 1) == 'g')
    {
        if (substring (ac_type, 3, 3) == 'a')
            ac_type <- 'geary'
        else
            ac_type <- 'getis_ord'
    } else
        ac_type <- 'moran'

    if (!missing (seed)) set.seed (seed)

    maxnbs <- max (sapply (nbs, length))
    nbsi <- lapply (seq (maxnbs), function (i) get_nbsi (i, nbs))

    if (!parallel)
    {
        dat <- rcpp_neutral_hotspots_ntests (nbs, wts, nbsi, alpha=alpha, sd0=sd0,
                                             niters=niters, ac_type=ac_type,
                                             log_scale=log_scale, ntests=ntests)
    } else
    {
        # Note that `makeCluster` with `type="FORK"` automatically attaches all
        # environment variables, but is not portable, as detailed at
        # r-bloggers.com/how-to-go-parallel-in-r-basics-tips/ 
        clust <- parallel::makeCluster (parallel::detectCores () - 1)
        exports <- list ('nbs', 'wts', 'nbsi', 'alpha', 'sd0',
                         'log_scale', 'niters', 'ac_type')
        parallel::clusterExport (clust, exports, envir=environment ())
        invisible (parallel::clusterCall (clust, function () {
                              while (length (grep ('hotspotr', getwd ())) > 0) 
                                  setwd ('..')
                              devtools::load_all ('hotspotr')
                              setwd ('./hotspotr')
                    }))

        z <- parallel::parLapply (clust, seq (ntests), function (i) 
                                  {
                                      rcpp_neutral_hotspots (nbs, wts, nbsi,
                                                             alpha=alpha,
                                                             sd0=sd0,
                                                             log_scale=log_scale, 
                                                             niters=niters,
                                                             ac_type=ac_type)
                                  })

        parallel::stopCluster (clust)
        ac <- colMeans (do.call (rbind, lapply (z, function (i) i [,2])))
        z <- colMeans (do.call (rbind, lapply (z, function (i) i [,1])))
        dat <- cbind (z, ac)
    }
    return (dat)
}

get_nbsi <- function (i, nbs)
{

    res <- lapply (seq (nbs), function (j)
                   {
                       if (length (nbs [[j]]) >= i)
                           c (j, nbs [[j]] [i], length (nbs [[j]]))
                       else
                           NULL
                   })
    res <- res [lapply (res, length) != 0]
    do.call (rbind, res)
}

#' order_one
#'
#' First order statistic for normal distribution.
#'
#' @param n Number of samples of normal distribution
#' @param sd Standard deviation of normal distribution
#' @param ntrials Number of trials over which to average order statistics
#'
#' @section Note Analytic calculation is possible following the first equation
#' from www.jstor.org/stable/2347982, as translated into R code adapted from
#' here: stackoverflow.com/questions/24211595/order-statistics-in-r
#'    integrand <- function (x, n, sigma=1) {
#'        x * pnorm (x, mean=1, sd=sigma, lower.tail=FALSE) ^ (n - 1) * 
#'            dnorm (x, mean=1, sd=sigma)
#'    }
#'
#'    o1 <- function(n, sigma=1) {
#'          integrate (integrand, -Inf, Inf, n, sigma)$value / beta (1, n)
#'    }
#' ... BUT the values of \code{o1} become wildly inaccurate for lower values of
#' \code{sd}, as can be seen, for example, for o1 (1e6, 0.01) = 7.2. This is
#' obviously nonsense, and the following numeric approximation is therefore
#' necessary.
#'
#' @return First order statistic
order_one <- function (n, sd, ntrials=1e4)
{
    if (missing (n)) stop ('n must be given')
    if (missing (sd)) stop ('sd must be given')

    temp <- rnorm (n=n, mean=1, sd=sd)
    if (min (temp) < 0)
        res <- 0
    else
        res <- mean (sapply (seq (ntrials), function (i) 
                             min (rnorm (n=n, mean=1, sd=sd))))

    return (res)
}
