#' neutral_hotspots2
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
neutral_hotspots2 <- function (nbs, wts, alpha = 0.1, sd0 = 0.1, ac_type = "moran",
                               niters = 1, log_scale = TRUE, ntests = 100, seed) {
    if (missing (nbs)) stop ("nbs must be given")

    if (missing (wts)) {
        wts <- lapply (nbs, function (x) rep (1, length (x)) / length (x))
    }

    # if (alpha [1] <= 0)
    #    stop ('neutral model only makes sense with finite temporal autocorrelation')

    ac_type <- tolower (ac_type)
    if (substring (ac_type, 1, 1) == "g") {
        if (substring (ac_type, 3, 3) == "a") {
            ac_type <- "geary"
        } else {
            ac_type <- "getis_ord"
        }
    } else {
        ac_type <- "moran"
    }

    if (!missing (seed)) set.seed (seed)

    maxnbs <- max (sapply (nbs, length))
    nbsi <- lapply (seq (maxnbs), function (i) get_nbsi (i, nbs))

    # number of sample points for density functions of environmental variables
    nvals <- 1e4
    x <- seq (0, 2, length.out = nvals)
    p0 <- truncnorm::dtruncnorm (x, a = 0, b = 2, mean = 1, sd = sd0)
    p1 <- array (p0, dim = c (nvals, length (nbs)))
    # expectation value:
    # integrate(f=function(x) truncnorm::dtruncnorm(x,a=0,b=2,mean=1,sd=0.5),
    #   lower=0,upper=2)
    # See
    # https://en.wikipedia.org/wiki/Product_distribution#Derivation
    # for the following calculations

    for (j in seq (niters))
    {
        p2 <- array (0, dim = c (nvals, length (nbs)))
        for (k in seq (maxnbs))
        {
            nbsi <- get_nbsi (k, nbs)
            nbsi <- data.frame (to = nbsi [, 1], from = nbsi [, 2], n = nbsi [, 3])
            p2 [, nbsi$to] <- p2 [, nbsi$to] + ((1 - alpha) * p1 [, nbsi$to] +
                alpha * p1 [, nbsi$from]) / nbsi$n
        }
        p1 <- p2
    }
    p1 <- rowMeans (p1)

    plot (x, p0, "l")
    lines (x, z1, col = "blue")

    if (log_scale) z1 <- log10 (z1)

    ac1 <- rcpp_ac_stats (z1, nbs, wts, ac_type)
    z1 <- (sort (z1, decreasing = TRUE) -
        min (z1)) / diff (range (z1))
    rbind (z1, ac1)


}
