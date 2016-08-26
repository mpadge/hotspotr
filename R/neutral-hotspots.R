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
#' @param niters Number of successive layers of spatial autocorrelation
#' @param log_scale If TRUE, raw hotspot values are log-transformed
#' @param ntests Number of tests over which to generate an average result
#' @param seed Random seed
#'
#' @return A vector of hotspot values sorted from high to low
#'
#' @examples
#' size <- 10
#' xy <- cbind (rep (seq (size), each=size), rep (seq (size), size))
#' dhi <- 1 # for rook; dhi=1.5 for queen
#' nbs <- spdep::dnearneigh (xy, 0, dhi)
#' z <- neutral_hotspots (nbs=nbs)
#'
#' @export
neutral_hotspots <- function (nbs, wts, alpha=0.1, sd0=0.1, niters=1, 
                              log_scale=TRUE, ntests=100, seed)
{
    if (missing (nbs)) stop ('nbs must be given')

    if (missing (wts)) 
        wts <- lapply (nbs, function (x) rep (1, length (x)) / length (x))

    if (alpha [1] <= 0)
        stop ('neutral model only makes sense with finite temporal autocorrelation')

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
        res <- do.call (rbind, res)
        data.frame (to=res [,1], from=res [,2], n=res [,3])
    }
    maxnbs <- max (sapply (nbs, length))

    z <- lapply (seq (ntests), function (i) 
                 {
                     z1 <- msm::rtnorm (size, mean=1, sd=sd0, lower=0, upper=2)
                     for (j in seq (niters))
                     {
                         z2 <- rep (0, size)
                         for (k in seq (maxnbs))
                         {
                             nbsi <- get_nbsi (k)
                             z2 [nbsi$to] <- z2 [nbsi$to] + 
                                 ((1 - alpha) * z1 [nbsi$to] +
                                  alpha * z1 [nbsi$from]) / nbsi$n
                         }
                         z1 <- z2
                     }
                     if (log_scale) z1 <- log10 (z1)
                     ac1 <- rcpp_ac_stats (z1, nbs, wts, ac_type)
                     z1 <- sort (z1, decreasing=TRUE)
                     z1 <- (z1 - min (z1)) / diff (range (z1))
                     cbind (z1, ac1)
                 })

    ac1 <- colMeans (do.call (rbind, lapply (z, function (i) i [,2])))
    z <- colMeans (do.call (rbind, lapply (z, function (i) i [,1])))
    data.frame (z=z, ac=ac1)
}
