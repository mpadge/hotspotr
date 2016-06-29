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
#' @param plot If TRUE, plot mean and observed distributions of z and associated
#' autocorrelation statistics
#' @param verbose If TRUE, dump progress details to screen
#'
#' @return Nothing (dumps statistics to screen)
#'
#' @export
p_values <- function (z, nbs, wts, alpha=c(0.1, 0.1), nt=100, ntests=1000,
                      ac_type='moran', plot=FALSE, verbose=FALSE)
{
    if (length (z) != length (nbs))
        stop ('z must be same size as nbs')

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

    # mean rank-scale distributions:
    if (verbose) 
        message ('Generating mean rank-scale distributions ... ',
                 appendLF=FALSE)
    rs_means <- neutral_hotspots_ntests (nbs=nbs, wts=wts, alpha=alpha, sd0=0.1,
                                         nt=nt, ntests=ntests, ac_type=ac_type)

    # Summed squared differences in rank-scale distributions:
    if (verbose) 
    {
        message ('done\nGenerating difference statistics for ', appendLF=FALSE)
        message ('rank-scale distributions ... ', appendLF=FALSE)
    }
    distributions <- rcpp_rs_dist_diff (nbs=nbs, wts=wts, alpha_t=alpha [1],
                                        alpha_s=alpha [2], sd0=0.1, nt=100,
                                        ntests=ntests, ac_type=ac_type,
                                        z_mn=rs_means [,1], ac_mn=rs_means [,2])
    if (verbose) message ('done')

    ac <- rcpp_ac_stats (z, nbs, wts, ac_type)
    z <- (sort (z, decreasing=TRUE) - min (z)) / diff (range (z))

    get_p <- function (z, zref, distrib)
    {
        stat <- sum ((z - zref) ^ 2)
        dens <- density (distrib)
        y <- cumsum (dens$y)
        p <- 1 - approx (dens$x, y / max (y), xout=stat)$y
        if (is.na (p)) p <- 0
        return (p)
    }
    p_z <- get_p (z, rs_means [,1], distributions [,1])
    p_ac <- get_p (ac, rs_means [,2], distributions [,2])

    if (plot)
    {
        plot.new ()
        par (mfrow=c(1,2))
        plot (1:length (z), z, "l", xlab="rank", ylab="scale")
        lines (1:length (z), rs_means [,1], col="gray")
        legend ("topright", lwd=1, col=c("black", "gray"), bty="n",
                legend=c("test data", "average distribution"))
        title (main=paste0 ("z: p = ", 
                            formatC (p_z, format="f", digits=4)))

        plot (1:length (ac), ac, "l", xlab="rank", ylab="scale")
        lines (1:length (ac), rs_means [,2], col="gray")
        title (main=paste0 ("ac: p = ", 
                            formatC (p_ac, format="f", digits=4)))
    }

    list (p_z=p_z, p_ac=p_ac)
}
