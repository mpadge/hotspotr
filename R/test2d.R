#' test2d
#'
#' Tests an observed matrix of data (\code{ymat}) against a two-dimensional
#' neutral model
#'
#' @param z Vector of observed values to be tested
#' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
#' point
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
#'   \item Temporal autocorrelation coefficient
#'   \item Spatial autocorrelation coefficient
#'   \item Number of successive layers of spatio-temporal autocorrelation
#'   required to reproduce statistical properties of observed data
#'   \item Absolute difference between observed and modelled rank--scale
#'   distributions
#' }
#'
#' @seealso \code{test1d}
#'
#' @examples
#' \dontrun{
#' alpha <- c (0.1, 0.1)
#' ymat <- ives2d (size=10, nt=1000, sd0=0.1, alpha=alpha)
#' result <- test2d (ymat, alpha=alpha, ntests=10)
#' }
#'
#' @export
test2d <- function (z, nbs, alpha=c(0.1, 0.1), ntests=100, ac_type='moran',
                    verbose=FALSE, plot=FALSE)
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

    ac <- rcpp_ac_stats (nbs, z, ac_type)
    zs <- sort ((z - min (z)) / diff (range (z)), decreasing=TRUE)
    test <- NULL # remove no visible binding warning
    # Initial 3D optimisation to get nt
    fn_n <- function (x)
    {
        test <- rcpp_neutral2d_ntests (nbs=nbs, alpha_t=x[1],
                                        alpha_s=x[2], sd0=0.1,
                                        nt=x[3], ntests=ntests, ac_type=ac_type)
        sum ((test [,1] - zs) ^ 2) + sum ((test [,2] - ac) ^ 2)
    }
    if (verbose) message ('Optimising for nt ... ', appendLF=FALSE)
    op <- optim (c (alpha, 10), fn_n)
    # then reduce to 2d optimisation
    nt <- round (op$par [3])
    alpha <- op$par [1:2]
    fn_a <- function (x)
    {
        test <- rcpp_neutral2d_ntests (nbs=nbs, alpha_t=x [1], alpha_s=x [2],
                                        sd0=0.1, nt=nt, ntests=ntests,
                                        ac_type=ac_type)
        sum ((test [,1] - zs) ^ 2) + sum ((test [,2] - ac) ^ 2)
    }
    if (verbose) message ('done.\nOptimising for alpha ... ', appendLF=FALSE)
    op <- optim (alpha, fn_a)
    if (verbose) message ('done.')
    alpha <- op$par
    test <- rcpp_neutral2d_ntests (nbs=nbs, alpha_t=alpha [1], alpha_s=alpha [2],
                                    sd0=0.1, nt=nt, ntests=ntests,
                                    ac_type=ac_type)

    # Parameters for ydat have been estimated; now generate equivalent neutral
    # values
    test_z <- as.numeric (test [,1])
    test_ac <- as.numeric (test [,2])

    # Note paired=TRUE is not appropriate because the positions in the sorted
    # lists are arbitrary and not directly related
    pval_z <- t.test (test_z, zs, paired=FALSE)$p.value
    pval_ac <- t.test (test_ac, ac, paired=FALSE)$p.value
    val_z <- sum ((test_z - zs) ^ 2)
    val_ac <- sum ((test_ac - ac) ^ 2)

    if (plot)
    {
        plot.new ()
        par (mfrow=c(1,2))
        cols <- c ('blue', 'red')
        y1 <- list (zs, ac)
        y2 <- list (test_z, test_ac)
        pvals <- c (pval_z, pval_ac)
        mt <- c ('raw', 'AC')
        for (i in 1:2)
        {
            plot (seq (y1 [[i]]), y1 [[i]], 'l', col=cols [1],
                  xlab='rank', ylab='scale')
            lines (seq (y2 [[i]]), y2 [[i]], col=cols [2])
            legend ('topright', lwd=1, col=cols, bty='n',
                    legend=c('observed', 'neutral2d'))
            title (main=paste0 (mt [i], ': p = ', 
                                formatC (pvals [i], format='f', digits=4)))
        }
    }

    pars <- list (alpha=alpha, nt=nt)
    pvals <- list (raw=pval_z, ac=pval_ac)
    data <- list (z=test_z, ac=test_ac)
    list (pars=pars, pvals=pvals, data=data)
}
