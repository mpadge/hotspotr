#' run_tests
#'
#' Tests observed data (\code{ydat}) against a series of neutral models
#'
#' @param size Size of the square grid on which to generate model. Total number
#' of points is size ^ 2. (Has no effect if \code{ydat} is given, in which case
#' value is determined by its dimensions)
#' @param alpha Vector of two components providing starting values for the
#' strength of autocorrelation in time and space
#' @param ntests Number of tests to run, with statistics calculated from the
#' mean of \code{ntests}
#' @param ydat Square matrix of observed values to be tested. If not given, will
#' be generated with specified values of \code{size} and \code{alpha}
#' @param seed Random seed
#'
#' @return Nothing (dumps statistics to screen)
#'
#' @seealso \code{test1d}, \code{test2d}
#'
#' @examples
#' \dontrun{
#' run_tests ()
#' }
#'
#' @export
run_tests <- function (size=10, alpha=c(0.1, 0.1), ntests=100, ydat, seed) 
{
    if (!missing (ydat))
    {
        if (!is.matrix (ydat)) stop ('ydat must be a matrix')
        if (dim (ydat) [1] != dim (ydat) [2]) stop ('ydat must be square')
    }

    # interesting seeds: 9
    if (missing (ydat))
    {
        if (missing (seed))
            ydat <- ives2D (size, 1000, sd0=0.1, alpha=alpha)
        else
            ydat <- ives2D (size, 1000, sd0=0.1, alpha=alpha,
                            seed=seed)
    }
    size <- sqrt (length (ydat))

    # To see how repeatable the tests are, they are performed with a different
    # seed
    set.seed (Sys.time ())
    
    cat ("  dim\t|\talpha\tdiff\tp(T)\t|\talpha\t\tn\t|\n", sep="")
    cat (rep ("-", 8), "|", rep ("-", 31), "|", rep ("-", 31), "|\n", sep="")
    alpha_s <- c (0.1, 0)
    for (i in 1:2)
    {
        t1 <- test1d (ydat, alpha=c(0.1, alpha_s [i]), ntests=ntests)
        cat ("  1\t|   (0.1, ", alpha_s [i], ")\t",
             formatC (sum ((t1$y - ydat) ^ 2), format="f", digits=2), "\t",
             formatC (t1$p.value, format="f", digits=4), "\t|   (",
             formatC (t1$pars$alpha [1], format="f", digits=2), ", ",
             formatC (t1$pars$alpha [2], format="f", digits=2), ")\t",
             round (t1$pars$n), "\t|\n", sep="")
    }

    alpha_t <- 0.1
    alpha_s <- c (0.1, 0)
    for (i in 1:2)
    {
        t2 <- test2d (ydat, alpha=c(alpha_t, alpha_s [i]), ntests=ntests)
        cat ("  2\t|   (0.1, ", alpha_s [i], ")\t",
             formatC (sum ((t2$y - ydat) ^ 2), format="f", digits=2), "\t",
             formatC (t2$p.value, format="f", digits=4), "\t|   (",
             formatC (t2$pars$alpha [1], format="f", digits=2), ", ",
             formatC (t2$pars$alpha [2], format="f", digits=2), ")\t",
             round (t2$pars$n), "\t|\n", sep="")
    }
    cat (rep ("-", 73), "\n", sep="")
}
