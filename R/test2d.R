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
test2d <- function (ydat, size=10, alpha=c(0.1, 0.1), sann=FALSE, plot=FALSE)
{
    ydat <- sort (ydat, decreasing=TRUE)
    ydat <- (ydat - min (ydat)) / diff (range (ydat))
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

    ytest <- neutral2d (size=size, alpha=op$par, n=n)
    ytest <- (ytest - min (ytest)) / diff (range (ytest))
    val <- sum ((ytest - ydat) ^ 2)

    if (plot)
    {
        plot (seq (ytest), ydat, "l", xlab="rank", ylab="scale")
        lines (seq (ytest), ytest, col="gray")
        legend ("topright", lwd=1, col=c("black", "gray"), bty="n",
                legend=c("observed", "neutral2d"))
    }

    c (op$par, n, val)
}
