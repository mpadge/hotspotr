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
#' ymat <- ives2D (size=10, nt=1000, sd0=0.1, alpha=alpha)
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
    fn_ac <- NULL
    if (substring (actype, 1, 1) == 'g')
    {
        if (substring (actype, 3, 3) == 'a')
            fn_ac <- 'gearyc'
        else
            fn_ac <- 'getis_ord'
    } else
        fn_ac <- 'morani'

    # Optimise until convergence based on raw values, with rank--scale
    # distributions for spatial AC statistics compared using the resultant
    # models. (Convergence for both raw data and AC statistics would be too
    # time-consuming.)
    size <- dim (ymat) [1]

    mdat <- do.call (fn_ac, list (ymat))
    ydat <- sort ((ymat - min (ymat)) / diff (range (ymat)), decreasing=TRUE)
    ytest <- NULL # remove no visible binding warning
    fn_n <- function (x)
    {
        ytest <- neutral2d (size=size, alpha=alpha, nt=x)
        ytest <- (ytest - min (ytest)) / diff (range (ytest))
        ytest <- sort (ytest, decreasing=TRUE)
        sum ((ytest - ydat) ^ 2)
    }
    nt <- round (optimise (fn_n, c (0, 200))$minimum)
    fn_a <- function (x)
    {
        ytest <- neutral2d (size=size, alpha=x, nt=nt)
        ytest <- (ytest - min (ytest)) / diff (range (ytest))
        ytest <- sort (ytest, decreasing=TRUE)
        sum ((ytest - ydat) ^ 2)
    }

    conv <- 1e99
    tol <- 1e-3
    nt0 <- c0 <- 100
    a0 <- a <- alpha
    count <- 0
    max_count <- 10
    while (conv > tol & count < max_count)
    {
        flags <- rep (FALSE, 2)
        op1 <- optimise (fn_n, c (1, 20))
        if (op1$objective < conv) 
        {
            conv <- op1$objective
            nt <- op1$minimum
            flags [1] <- TRUE
        }
        op <- optim (alpha, fn_a)
        if (op$value < conv) 
        {
            conv <- op$value
            alpha <- op$par
            flags [2] <- TRUE
        }

        if (!all (flags) & c0 == conv) 
            count <- count + 1
        else
            count <- 0
        a0 <- alpha
        nt0 <- nt
        c0 <- conv
    }

    # Parameters for ydat have been estimated; now generate equivalent neutral
    # values
    yt <- ym <- rep (0, size ^ 2)
    for (i in 1:ntests)
    {
        yi <- neutral2d (size, alpha=a0, nt=nt0)
        yt <- yt + sort ((yi - min (yi)) / diff (range (yi)), decreasing=TRUE)
        mi <- do.call (fn_ac, list (yi))
        ym <- ym + sort ((mi - min (mi)) / diff (range (mi)), decreasing=TRUE)
    }
    yt <- yt / ntests
    ym <- ym / ntests
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
