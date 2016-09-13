#' generate_hotspot_model
#'
#' Uses the parameters returned by \code{fit_hotspot_model} to generate
#' rank-scale distributions both of raw values and associated spatial
#' autocorrelation statistics.
#'
#' @param n Number of observations to generate
#' @param alpha Strength of spatial autocorrelation 
#' @param sd0 Standard deviation of truncated normal distribution used to model
#' environmental variation (with mean of 1)
#' @param ac_type Type of autocorrelation statistic to use in tests
#' (\code{moran}, \code{geary}, or \code{getis-org}=\code{go})
#' @param niters Number of successive layers of spatial autocorrelation
#' @param plot If TRUE, produces a plot of rank--scale distributions
#'
#' @return A matrix of two columns containing sorted and scaled versions of
#' \enumerate{
#'   \item z = raw values
#'   \item ac = associated spatial autocorrelation statistics
#' }
#'
#' @export
generate_hotspot_model <- function (n, alpha=0.1, sd0=0.1, ac_type='moran',
                                    niters=1, plot=FALSE)
{
    z1 <- msm::rtnorm (n, mean=1, sd=sd0, lower=0, upper=2)
    for (j in seq (niters))
    {
        z2 <- rep (0, size)
        for (k in seq (maxnbs))
        {
            nbsi <- get_nbsi (k)
            z2 [nbsi$to] <- z2 [nbsi$to] + 
                ((1 - alpha) * z1 [nbsi$to] +
                 alpha * z1 [nbsi$from]) / nbsi$n
        }
        z1 <- z2
    }
    if (log_scale) z1 <- log10 (z1)
    ac1 <- rcpp_ac_stats (z1, nbs, wts, ac_type)
    z1 <- sort (z1, decreasing=TRUE)
    z1 <- (z1 - min (z1)) / diff (range (z1))
    cbind (z1, ac1)

}
