#' test1d
#'
#' Tests an observed matrix of data (\code{ydat}) against a one-dimensional
#' neutral model
#'
#' @param ydat Square matrix of observed values to be tested
#' @param alpha Vector of two components providing starting values for the
#' strength of autocorrelation in time and space
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
#' @seealso \code{test2d}
#'
#' @export
test1d <- function (ydat, alpha=c(0.1, 0.1))
{
    # progressive optimisation until convergence
    size <- sqrt (length (ydat))

    ydat <- sort (ydat, decreasing=TRUE)
    ydat <- (ydat - min (ydat)) / diff (range (ydat))
    fn_n <- function (x)
    {
        ytest <- neutral1d (size, alpha=alpha, n=x)
        ytest <- (ytest - min (ytest)) / diff (range (ytest))
        sum ((ytest - ydat) ^ 2)
    }
    fn_a <- function (x)
    {
        ytest <- neutral1d (size, alpha=x, n=n)
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
        op1 <- optimise (fn_n, c (0, 200))
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
    ytest <- neutral1d (size, alpha=alpha, n=n)
    ytest <- (ytest - min (ytest)) / diff (range (ytest))
    val <- sum ((ytest - ydat) ^ 2)

    c (alpha, n, val)
}
