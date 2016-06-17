#' test2d
#'
#' Tests an observed matrix of data (\code{ydat}) against a two-dimensional
#' neutral model
#'
#' @param ydat Square matrix of observed values to be tested
#' @param alpha Vector of two components providing starting values for the
#' strength of autocorrelation in time and space
#' @param ntests Number of repeats of neutral model used to calculate mean
#' rank--scale distribution
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
#' ydat <- ives2D (size=10, nt=1000, sd0=0.1, alpha=alpha)
#' result <- test2d (ydat, alpha=alpha, ntests=10)
#' }
#'
#' @export
test2d <- function (ydat, alpha=c(0.1, 0.1), ntests=100, 
                    plot=FALSE)
{
    size <- sqrt (length (ydat))
    ydat <- sort (ydat, decreasing=TRUE)
    ydat <- (ydat - min (ydat)) / diff (range (ydat))
    ytest <- NULL # remove no visible binding warning
    fn_n <- function (x)
    {
        ytest <- neutral2d (size=size, alpha=alpha, n=x)
        ytest <- (ytest - min (ytest)) / diff (range (ytest))
        sum ((ytest - ydat) ^ 2)
    }
    n <- round (optimise (fn_n, c (0, 200))$minimum)
    fn_a <- function (x)
    {
        ytest <- neutral2d (size=size, alpha=x, n=n)
        ytest <- (ytest - min (ytest)) / diff (range (ytest))
        sum ((ytest - ydat) ^ 2)
    }

    conv <- 1e99
    tol <- 1e-3
    n0 <- c0 <- 100
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
            n <- op1$minimum
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
        n0 <- n
        c0 <- conv
    }

    # Parameters for ydat have been estimated; now generate equivalent neutral
    # values
    ytest <- rep (0, size ^ 2)
    for (i in 1:ntests)
        ytest <- ytest + neutral2d (size, alpha=a0, n=n0)
    ytest <- ytest / ntests
    # Note paired=TRUE is not appropriate because the positions in the sorted
    # lists are arbitrary and not directly related
    pval <- t.test (ytest, ydat, paired=FALSE)$p.value
    val <- sum ((ytest - ydat) ^ 2)

    if (plot)
    {
        cols <- c ('blue', 'red')
        plot (seq (ytest), ydat, 'l', xlab='rank', ylab='scale', col=cols [1])
        lines (seq (ytest), ytest, col=cols [2])
        legend ('topright', lwd=1, col=cols, bty='n',
                legend=c('observed', 'neutral2d'))
        title (main=paste0 ('p = ', formatC (pval, format='f', digits=4)))
    }

    pars <- list (alpha=a0, n=n)
    list (pars=pars, difference=val, p.value=pval, y=ytest)
}
