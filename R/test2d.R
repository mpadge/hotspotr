#' test2d
#'
#' Tests an observed matrix of data (\code{ydat}) against a two-dimensional
#' neutral model
#'
#' @param ydat Square matrix of observed values to be tested
#' @param size Size of the square grid on which to generate model. Total number
#' of points is size ^ 2
#' @param alpha Vector of two components providing starting values for the
#' strength of autocorrelation in time and space
#' @param ntests Number of repeats of neutral model used to calculate mean
#' rank--scale distribution
#' @param sann Use simulated annealing to find optimum
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
#' @export
test2d <- function (ydat, size=10, alpha=c(0.1, 0.1), ntests=100, sann=FALSE, plot=FALSE)
{
    ydat <- sort (ydat, decreasing=TRUE)
    ydat <- (ydat - min (ydat)) / diff (range (ydat))
    ytest <- NULL # remove no visible binding warning
    fn1 <- function (x)
    {
        ytest <- neutral2d (size=size, alpha=alpha, n=x)
        ytest <- (ytest - min (ytest)) / diff (range (ytest))
        sum ((ytest - ydat) ^ 2)
    }
    n <- round (optimise (fn1, c (0, 200))$minimum)
    fn2 <- function (x)
    {
        ytest <- neutral2d (size=size, alpha=x, n=n)
        ytest <- (ytest - min (ytest)) / diff (range (ytest))
        sum ((ytest - ydat) ^ 2)
    }

    if (sann)
        op <- optim (alpha, fn2, method='SANN')
    else
        op <- optim (alpha, fn2)

    # Parameters for ydat have been estimated; now generate equivalent neutral
    # values
    dd <- rep (NA, ntests)
    y1s <- rep (0, size ^ 2)
    for (j in 1:ntests)
    {
        y1 <- neutral1d (size=size, alpha=op$par, n=n)
        y1 <- (y1 - min (y1)) / diff (range (y1))
        dd [j] <- sum (y1 - ydat) 
        y1s <- y1s + y1
    }
    y1s <- y1s / ntests
    pval <- t.test (dd)$p.value
    val <- sum ((y1s - ydat) ^ 2)

    if (plot)
    {
        plot (seq (ytest), ydat, "l", xlab="rank", ylab="scale")
        lines (seq (ytest), ytest, col="gray")
        legend ("topright", lwd=1, col=c("black", "gray"), bty="n",
                legend=c("observed", "neutral2d"))
    }

    pars <- list (alpha=op$par, n=n)
    list (pars=pars, difference=val, p.value=pval, y=y1s)
}
