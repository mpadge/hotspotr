#' fit_hotspot_model
#'
#' Fits a neutral model to an observed series of hotspot values in terms of
#' temporal and spatial autocorrelation parameters, and number of iterations of
#' these.
#'
#' @param z Vector of observed values to be tested
#' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
#' point
#' @param wts Weighting factors for each neighbour; must have same length as
#' nbs. Uniform weights used if not given.
#' @param alpha Vector of two components providing starting values for the
#' strength of autocorrelation in time and space
#' @param ntests Number of repeats of neutral model used to calculate mean
#' rank--scale distribution
#' @param ac_type type of autocorrelation statistic to use in tests
#' (\code{moran}, \code{geary}, or \code{getis-ord}=\code{go})
#' @param verbose If TRUE, dump progress details to screen
#' @param plot If TRUE, produces a plot of rank--scale distributions
#'
#' @return A vector of four values as estimated by the neutral model:
#' \enumerate{
#'   \item alpha respectively containing temporal and spatial autocorrelation
#'   coefficients
#'   \item nt = Number of successive layers of spatio-temporal autocorrelation
#'   required to reproduce statistical properties of observed data
#' }
#'
#' @examples
#' \dontrun{
#' alpha <- c (0.1, 0.1)
#' dat <- ives (size=10, nt=1000, sd0=0.1, alpha=alpha)
#' test <- fit_hotspot_model (nbs=dat$nbs, z=dat$z, alpha=alpha, ntests=10)
#' }
#'
#' @export
fit_hotspot_model <- function (z, nbs, wts, alpha=c(0.1, 0.1), ntests=100,
                               ac_type='moran', verbose=FALSE, plot=FALSE)
{
    if (!is.numeric (z)) 
        stop ('z must be numeric')
    if (!is (nbs, 'nb'))
        stop ('nbs must of class spdep::nb')
    if (length (z) != length (nbs))
        stop ('nbs must have same length as z')

    ac_type <- tolower (ac_type)
    if (substring (ac_type, 1, 1) == 'g')
    {
        if (substring (ac_type, 3, 3) == 'a')
            ac_type <- 'geary'
        else
            ac_type <- 'getis_ord'
    } else
        ac_type <- 'moran'

    size <- length (z)

    if (missing (wts)) 
        wts <- lapply (nbs, function (x) rep (1, length (x)) / length (x))

    ac <- rcpp_ac_stats (z, nbs, wts, ac_type)
    zs <- sort ((z - min (z)) / diff (range (z)), decreasing=TRUE)
    test <- NULL # remove no visible binding warning
    # Initial 3D optimisation to get nt
    fn_n <- function (x)
    {
        test <- rcpp_neutral_hotspots_ntests (nbs=nbs, wts=wts, alpha_t=x[1],
                                              alpha_s=x[2], sd0=0.1,
                                              nt=x[3], ntests=ntests, 
                                              ac_type=ac_type)
        sum ((test [,1] - zs) ^ 2) + sum ((test [,2] - ac) ^ 2)
    }
    control <- list (reltol=1e-3)
    if (verbose) 
    {
        message ('Optimising for nt ... ', appendLF=FALSE)
        #control <- c (control, trace=1)
    }
    op <- optim (c (alpha, 10), fn_n, control=control)
    # then reduce to 2d optimisation
    nt <- round (op$par [3])
    alpha <- op$par [1:2]
    fn_a <- function (x)
    {
        test <- rcpp_neutral_hotspots_ntests (nbs=nbs, wts=wts, alpha_t=x [1],
                                              alpha_s=x [2], sd0=0.1, nt=nt,
                                              ntests=ntests, ac_type=ac_type)
        sum ((test [,1] - zs) ^ 2) + sum ((test [,2] - ac) ^ 2)
    }
    if (verbose) message ('done.\nOptimising for alpha ... ', appendLF=FALSE)
    op <- optim (alpha, fn_a, control=control)
    if (verbose) message ('done.')
    alpha <- op$par

    list (alpha=alpha, nt=nt)
}
