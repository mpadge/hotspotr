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
                                     ac_type='moran', seed)
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

    # single vector indices of all neighbours of each points and associated wts
    indx_nbso <- unlist (lapply (seq (nbs), function (i) 
                                 rep (i, length (nbs [[i]]))))
    indx_nbs <- unlist (nbs)
    indx_wts <- unlist (wts)

    size <- length (nbs)

    # Set up parallel cluster:
    clust <- parallel::makeCluster (parallel::detectCores () - 1)
    # "envir" necessary to get variables from this environment
    #parallel::clusterExport (clust, list ("nbs", "wts", "alpha", "sd0", "nt",
    #                                      "ac_type"), envir=environment ())
    exports <- list ("size", "nbs", "wts", "ac_type", "alpha", "sd0",
                     "ntests", "indx_nbso", "indx_nbs", "indx_wts")
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
                                  z1 [indx_nbso] <- (1 - alpha) * 
                                      z1 [indx_nbso] + alpha * z1 [indx_nbs] * 
                                      indx_wts [indx_nbs]
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
