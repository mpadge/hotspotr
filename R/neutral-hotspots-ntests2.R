#' neutral_hotspots_ntests2
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
#' \dontrun{
#' nbs <- ives (size=10)$nbs
#' z <- neutral_hotspots (nbs=nbs)
#' }
#'
#' @export
neutral_hotspots_ntests2 <- function (nbs, wts, alpha=0.1, sd0=0.1, niters=1, 
                                      log_scale=TRUE, ntests=100, parallel=FALSE,
                                      seed)
{
    if (missing (nbs)) stop ('nbs must be given')

    if (missing (wts)) 
        wts <- lapply (nbs, function (x) rep (1, length (x)) / length (x))

    #if (alpha [1] <= 0)
    #    stop ('neutral model only makes sense with finite temporal autocorrelation')

    if (!missing (seed)) set.seed (seed)

    size <- length (nbs)

    ac_type <- 'moran'

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
        do.call (rbind, res)
    }
    maxnbs <- max (sapply (nbs, length))
    nbsi <- lapply (seq (maxnbs), function (i) get_nbsi (i))

    if (!parallel)
    {
        rcpp_neutral_hotspots_ntests (nbs, wts, nbsi, alpha=alpha, sd0=sd0,
                                      niters=niters, ac_type=ac_type,
                                      log_scale=log_scale, ntests=ntests)
    } else
    {
        # Note that `makeCluster` with `type="FORK"` automatically attaches all
        # environment variables, but is not portable, as detailed at
        # r-bloggers.com/how-to-go-parallel-in-r-basics-tips/ 
        clust <- parallel::makeCluster (parallel::detectCores () - 1)
        exports <- list ("nbs", "wts", "nbsi", "alpha", "sd0",
                         "log_scale", "niters", "ac_type")
        parallel::clusterExport (clust, exports, envir=environment ())
        invisible (parallel::clusterCall (clust, function () {
                              while (length (grep ('hotspotr', getwd ())) > 0) 
                                  setwd ("..")
                              devtools::load_all ("hotspotr")
                              setwd ("./hotspotr")
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
        cbind (z, ac)
    }
}
