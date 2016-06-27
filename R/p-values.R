#' p_values
#'
#' Tests observed data (\code{z}) against a series of neutral models
#'
#' @param z Vector of observed values to be tested
#' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
#' point
#' @param wts Weighting factors for each neighbour; must have same length as
#' nbs. Uniform weights used if not given.
#' @param alpha Vector of two components providing starting values for the
#' strength of autocorrelation in time and space
#' @param nt Number of successive layers of temporal and spatial autocorrelation
#' used to generate test field
#' @param ntests Number of tests to run, with statistics calculated from the
#' mean of \code{ntests}
#' @param ac_type type of autocorrelation statistic to use in tests
#' (\code{moran}, \code{geary}, or \code{getis-ord}=\code{go})
#'
#' @return Nothing (dumps statistics to screen)
#'
#' @export
p_values <- function (z, nbs, wts, alpha=c(0.1, 0.1), nt=100, ntests=1000,
                      ac_type='moran')
{
    if (length (z) != length (nbs))
        stop ('z must be same size as nbs')

    ac_type <- tolower (ac_type)
    if (substring (ac_type, 1, 1) == 'g')
    {
        if (substring (ac_type, 3, 3) == 'a')
            ac_type <- 'geary'
        else
            ac_type <- 'getis_ord'
    } else
        ac_type <- 'moran'

    distributions <- rcpp_p_values (nbs=nbs, wts=wts, alpha_t=alpha [1],
                                    alpha_s=alpha [2], sd0=0.1, nt=100,
                                    ntests=ntests, ac_type=ac_type)
    # The distribution is **close to** but not the same as a chi-squared:
    # plot.new ()
    # par (mfrow=c(1,2))
    # mts <- c (expression (paste (z^2)), expression (paste (ac^2)))
    # mults <- c (100, 1000)
    # for (i in 1:2)
    # {
    #     y <- sqrt (test [,i]) # SDs
    #     x <- seq (0, ceiling (max (y) * mults [i]) / mults [i], by=1/mults [i])
    #     hh <- hist (y, breaks=x, plot=FALSE)
    #     plot (x [1:length (x) - 1] * mults [i], hh$counts, "l", 
    #           xlab="x", ylab="freq", main=mts [i])
    #
    #     x <- x [1:length (x) - 1] * mults [i]
    #     y <- dchisq (x, 21) # arbitrary value that fits observations well
    #     ry <- range (y, finite=TRUE)
    #     y <- (y - ry [1]) / diff (ry) * max (hh$counts)
    #     lines (x, y, col="red")
    #}
    # mean rank-scale distributions:
    rs_means <- rcpp_neutral_hotspots_ntests (nbs=nbs, wts=wts, alpha_t=alpha [1],
                                          alpha_s=alpha [2], sd0=0.1,
                                          nt=100, ntests=100, 
                                          ac_type=ac_type)
    z <- (sort (z, decreasing=TRUE) - min (z)) / diff (range (z))
    stat_z <- sum ((z - rs_means [,1]) ^ 2)
    p_z <- length (which (distributions [,1] > stat_z)) / ntests
    ac <- rcpp_ac_stats (z, nbs, wts, ac_type)
    stat_ac <- sum ((ac - rs_means [,2]) ^ 2)
    p_ac <- length (which (distributions [,2] > stat_ac)) / ntests

    list (p_z=p_z, p_ac=p_ac)
}
