#' test2d
#'
#' Tests an observed matrix of data (\code{ymat}) against a two-dimensional
#' neutral model
#'
#' @param ymat Square matrix of observed values to be tested
#' @param alpha Vector of two components providing starting values for the
#' strength of autocorrelation in time and space
#' @param ntests Number of repeats of neutral model used to calculate mean
#' rank--scale distribution
#' @param actype type of autocorrelation statistic to use in tests
#' (\code{moran}, \code{geary}, or \code{getis-ord}=\code{go})
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
test2d <- function (ymat, alpha=c(0.1, 0.1), ntests=100, actype='moran',
                    plot=FALSE)
{
    if (!is.matrix (ymat)) 
        stop ('ymat must be a matrix')
    else if (dim (ymat) [1] != dim (ymat) [2])
        stop ('ymat must be a square matrix')

    actype <- tolower (actype)
    fn_ac <- 'morani'
    ac_num <- 1
    if (substring (actype, 1, 1) == 'g')
    {
        if (substring (actype, 3, 3) == 'a')
        {
            fn_ac <- 'gearyc'
            ac_num <- 2
        } else
        {
            fn_ac <- 'getis_ord'
            ac_num <- 3
        }
    } 

    # Optimise using raw values only
    size <- dim (ymat) [1]

    mdat <- do.call (fn_ac, list (ymat))
    ydat <- sort ((ymat - min (ymat)) / diff (range (ymat)), decreasing=TRUE)
    yt <- NULL # remove no visible binding warning
    # Initial 3D optimisation to get nt
    fn_n <- function (x)
    {
        yt <- rcpp_neutral2d_ntests (size=size, alpha_t=x[1],
                                        alpha_s=x[2], sd0=0.1,
                                        nt=x[3], ntests=ntests, ac_type=ac_num)
        sum ((yt [,1] - ydat) ^ 2) # [,1] = raw values; [,2] = AC stats
    }
    op <- optim (c (alpha, 10), fn_n)
    # then reduce to 2d optimisation
    nt <- round (op$par [3])
    alpha <- op$par [1:2]
    fn_a <- function (x)
    {
        yt <- rcpp_neutral2d_ntests (size=size, alpha_t=x [1], alpha_s=x [2],
                                        sd0=0.1, nt=nt, ntests=ntests,
                                        ac_type=ac_num)
        sum ((yt [,1] - ydat) ^ 2)
    }
    op <- optim (alpha, fn_a)
    alpha <- op$par
    yt <- rcpp_neutral2d_ntests (size=size, alpha_t=alpha [1], alpha_s=alpha [2],
                                    sd0=0.1, nt=nt, ntests=ntests,
                                    ac_type=ac_num)

    # Parameters for ydat have been estimated; now generate equivalent neutral
    # values
    ym <- as.numeric (yt [,2])
    yt <- as.numeric (yt [,1])

    # Note paired=TRUE is not appropriate because the positions in the sorted
    # lists are arbitrary and not directly related
    pval_t <- t.test (yt, ydat, paired=FALSE)$p.value
    pval_m <- t.test (ym, mdat, paired=FALSE)$p.value
    val_t <- sum ((yt - ydat) ^ 2)
    val_m <- sum ((ym - ydat) ^ 2)

    if (plot)
    {
        plot.new ()
        par (mfrow=c(1,2))
        cols <- c ('blue', 'red')
        y1 <- list (ydat, mdat)
        y2 <- list (yt, ym)
        pvals <- c (pval_t, pval_m)
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
    pvals <- list (raw=pval_t, ac=pval_m)
    data <- list (y=yt, ac=ym)
    list (pars=pars, pvals=pvals, data=data)
}
