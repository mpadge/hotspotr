#' rs_dist_diff
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
#' @param mean_stats Mean rank--scale distributions returned from
#' \code{neutral_hotspots_ntests}: a matrix of 2 columns
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
rs_dist_diff <- function (nbs, wts, alpha=c(0.1, 0.1), sd0=0.1, ntests=1000,
                                     nt=100, ac_type='moran', mean_stats)
{
    if (missing (nbs)) stop ('nbs must be given')

    if (missing (wts)) 
        wts <- lapply (nbs, function (x) rep (1, length (x)) / length (x))

    if (missing (mean_stats)) 
        mean_stats <- neutral_hotspots_ntests (nbs, wts, alpha, ntests, nt,
                                               sd0, ac_type)
    z_mn <- mean_stats$z
    ac_mn <- mean_stats$ac

    # Set up parallel cluster:
    clust <- parallel::makeCluster (parallel::detectCores () - 1)
    # "envir" necessary to get variables from this environment
    parallel::clusterExport (clust, list ("nbs", "wts", "alpha", "sd0", "nt",
                                          "z_mn", "ac_mn", "ac_type"),
                             envir=environment ())
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
                        c (sum ((z1 - z_mn) ^ 2), sum ((ac1 - ac_mn) ^ 2))
                    })
    parallel::stopCluster (clust)
    z <- do.call (rbind, z)
    data.frame (z=z [,1], ac=z [,2])
}
