#' test1d
#'
#' Tests an observed matrix of data (\code{ydat}) against a one-dimensional
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
#' @seealso \code{test2d}
#'
#' @examples
#' \dontrun{
#' alpha <- c (0.1, 0.1)
#' ydat <- ives2D (size=10, nt=1000, sd0=0.1, alpha=alpha)
#' result <- test1d (ydat, alpha, ntests=10)
#' }
#'
#' @export
test1d <- function (ydat, alpha=c(0.1, 0.1), ntests=100, plot=FALSE)
{
    # progressive optimisation until convergence
    size <- sqrt (length (ydat))

    ydat <- sort (ydat, decreasing=TRUE)
    ydat <- (ydat - min (ydat)) / diff (range (ydat))
    ytest <- NULL # remove no visible binding warning
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
    # values. The whole routine with the following 2 lapply lines took 1.55s,
    # while the for-loop version took 1.61s.
    #ytest <- lapply (seq (ntests), function (i) neutral1d (size, alpha=a0, n=n0))
    #ytest <- colSums (do.call (rbind, ytest)) / ntests
    ytest <- rep (0, size ^ 2)
    for (i in 1:ntests)
        ytest <- ytest + neutral1d (size, alpha=a0, n=n0)
    ytest <- ytest / ntests
    # Note paired=TRUE is not appropriate because the positions in the sorted
    # lists are arbitrary and not directly related
    pval <- t.test (ytest, ydat, paired=FALSE)$p.value
    val <- sum ((ytest - ydat) ^ 2)

    if (plot)
    {
        plot (seq (ytest), ydat, "l", xlab="rank", ylab="scale")
        lines (seq (ytest), ytest, col="gray")
        legend ("topright", lwd=1, col=c("black", "gray"), bty="n",
                legend=c("observed", "neutral1d"))
        title (main=paste0 ("p = ", formatC (pval, format="f", digits=4)))
    }

    pars <- list (alpha=alpha, n=n)
    list (pars=pars, difference=val, p.value=pval, y=ytest)
}
