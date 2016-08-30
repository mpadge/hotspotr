#' rs_dist_diff
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
#' @param log_scale If TRUE, raw hotspot values are log-transformed
#' @param niters Number of sequential iterations of spatial autocorrelation
#' @param ntests Number of repeats of neutral model used to calculate mean
#' rank--scale distribution
#' @param ac_type type of autocorrelation statistic to use in tests
#' (\code{moran}, \code{geary}, or \code{getis-ord}=\code{go})
#' @param mean_stats Mean rank--scale distributions returned from
#' \code{neutral_hotspots}: a matrix of 2 columns
#'
#' @return A vector of hotspot values sorted from high to low
#'
#' @seealso \code{ives}
#'
#' @export
rs_dist_diff <- function (nbs, wts, alpha=0.1, sd0=0.1, niters=1, ntests=1000,
                          ac_type='moran', log_scale=TRUE, mean_stats)
{
    if (missing (nbs)) stop ('nbs must be given')

    if (missing (wts)) 
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

    maxnbs <- max (sapply (nbs, length))
    nbsi <- lapply (seq (maxnbs), function (i) get_nbsi (i, nbs))

    nt <- NULL # suppress no visisble binding package check note
    if (missing (mean_stats)) 
        mean_stats <- neutral_hotspots (nbs, wts, alpha=alpha, sd0=sd0,
                                        niters=niters, ac_type=ac_type,
                                        ntests=ntests)
    z_mn <- mean_stats [,1]
    ac_mn <- mean_stats [,2]

    # Individual tests are run across a parallel cluster with code exactly as in
    # neutral-hotspots.R
    clust <- parallel::makeCluster (parallel::detectCores () - 1)
    exports <- list ('nbs', 'wts', 'nbsi', 'alpha', 'sd0', 'log_scale',
                     'niters', 'ac_type', 'z_mn', 'ac_mn')
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
                                                         ac_type=ac_type) -
                                  cbind (z_mn, ac_mn)
                              })
    parallel::stopCluster (clust)
    ac <- do.call (rbind, lapply (z, function (i) i [,2]))
    z <- do.call (rbind, lapply (z, function (i) i [,1]))

    # result is two columns of summed squared differences for
    cbind (rowSums (z ^ 2), rowSums (ac ^ 2))
}
