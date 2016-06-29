#' neutral_hotspots_ntests
#'
#' Averages results of a number of neutral models
#'
#' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
#' point
#' @param wts Weighting factors for each neighbour; must have same length as
#' nbs. Uniform weights used if not given.
#' @param alpha Vector of two components respectively specifying the strength of
#' autocorrelation in time and space.
#' @param ntests Number of repeats of neutral model used to calculate mean
#' rank--scale distribution
#' @param nt Number of successive layers of temporal and spatial autocorrelation
#' used to generate final modelled values
#' @param sd0 Standard deviation of truncated normal distribution used to model
#' environmental variation (with mean of 1)
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
neutral_hotspots_ntests <- function (nbs, wts, alpha=c(0.1, 0.1), ntests=1000,
                                     nt=100, sd0=0.1, ac_type='moran', seed)
{
    if (missing (nbs)) stop ('nbs must be given')

    if (missing (wts)) 
        wts <- lapply (nbs, function (x) rep (1, length (x)) / length (x))

    if (!missing (seed)) set.seed (seed)

    # Set up parallel cluster:
    clust <- parallel::makeCluster (parallel::detectCores () - 1)
    parallel::clusterExport (clust, "nbs")
    parallel::clusterExport (clust, "wts")
    parallel::clusterExport (clust, "ac_type")
    invisible (parallel::clusterCall (clust, function () {
                            while (length (grep ('hotspotr', getwd ())) > 0) 
                                setwd ("..")
                            devtools::load_all ("hotspotr")
                            setwd ("./hotspotr")
                                             }))

    z <- parallel::parLapply (clust, seq (ntests), function (i) 
                    {
                        z1 <- rcpp_neutral_hotspots (nbs, wts, 
                                                     alpha_t=alpha [1], 
                                                     alpha_s=alpha [2], sd0=sd0,
                                                     nt=nt)
                        ac1 <- rcpp_ac_stats (z1, nbs, wts, ac_type)
                        z1 <- (sort (z1, decreasing=TRUE) - min (z1)) / 
                                diff (range (z1)) 
                        rbind (z1, ac1)
                    })
    parallel::stopCluster (clust)
    ac <- colMeans (do.call (rbind, lapply (z, function (i) i [2,])))
    z <- colMeans (do.call (rbind, lapply (z, function (i) i [1,])))
    data.frame (z=z, ac=ac)
}
