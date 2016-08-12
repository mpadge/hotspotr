#' neutral_hotspots_ntests
#'
#' Averages results of a number of neutral models
#'
#' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
#' point
#' @param wts Weighting factors for each neighbour; must have same length as
#' nbs. Uniform weights used if not given.
#' @param alpha Strength of spatial autocorrelation 
#' @param sd0 Standard deviation of truncated normal distribution used to model
#' environmental variation (with mean of 1)
#' @param niters Number of sequential iterations of spatial autocorrelation
#' @param ntests Number of repeats of neutral model used to calculate mean
#' rank--scale distribution
#' @param ac_type type of autocorrelation statistic to use in tests
#' (\code{moran}, \code{geary}, or \code{getis-ord}=\code{go})
#' @param seed Random seed
#'
#' @return A vector of hotspot values sorted from high to low
#'
#' @seealso \code{ives}
#'
#' @examples
#' nbs <- ives (size=10)$nbs
#' z <- neutral_hotspots (nbs=nbs)
#'
#' @export
neutral_hotspots_ntests <- function (nbs, wts, alpha=0.1, sd0=0.1, ntests=1000,
                                     niters=1, ac_type='moran', seed)
{
    if (missing (nbs)) stop ('nbs must be provided')
    if (!is (nbs, 'nb')) stop ('nbs must of class spdep::nb')
    if (!missing (wts))
    {
        if (length (wts) != length (nbs))
            stop ('wts must have the same length as nbs')
    } else
        wts <- lapply (nbs, function (x) rep (1, length (x)) / length (x))

    if (!missing (seed)) set.seed (seed)

    size <- length (nbs)

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

    # Set up parallel cluster:
    clust <- parallel::makeCluster (parallel::detectCores () - 1)
    # "envir" necessary to get variables from this environment
    #parallel::clusterExport (clust, list ("nbs", "wts", "alpha", "sd0", "nt",
    #                                      "ac_type"), envir=environment ())
    exports <- list ("size", "nbs", "wts", "ac_type", "alpha", "sd0",
                     "ntests", "get_nbsi", "maxnbs")
    parallel::clusterExport (clust, exports, envir=environment ())
    invisible (parallel::clusterCall (clust, function () {
                            while (length (grep ('hotspotr', getwd ())) > 0) 
                                setwd ("..")
                            devtools::load_all ("hotspotr")
                            setwd ("./hotspotr")
                                             }))

    #z <- parallel::parLapply (clust, seq (ntests), function (i) 
    #                {
    #                    z1 <- rcpp_neutral_hotspots (nbs, wts, 
    #                                                 alpha_t=alpha [1], 
    #                                                 alpha_s=alpha [2], sd0=sd0,
    #                                                 nt=nt)
    #                    ac1 <- rcpp_ac_stats (z1, nbs, wts, ac_type)
    #                    z1 <- (sort (z1, decreasing=TRUE) - min (z1)) / 
    #                            diff (range (z1)) 
    #                    rbind (z1, ac1)
    #                })

    z <- parallel::parLapply (clust, seq (ntests), function (i) 
                              {
                                  z1 <- msm::rtnorm (size, mean=1, sd=sd0,
                                                     lower=0, upper=2)

                                  for (j in seq (niters))
                                  {
                                      z2 <- rep (0, size ^ 2)
                                      for (k in seq (maxnbs))
                                      {
                                          nbsi <- get_nbsi (k)
                                          z2 [nbsi$to] <- z2 [nbsi$to] + 
                                              ((1 - alpha) * z1 [nbsi$to] +
                                               alpha * z1 [nbsi$from]) / nbsi$n
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

    parallel::stopCluster (clust)
    ac <- colMeans (do.call (rbind, lapply (z, function (i) i [2,])))
    z <- colMeans (do.call (rbind, lapply (z, function (i) i [1,])))
    data.frame (z=z, ac=ac)
}
