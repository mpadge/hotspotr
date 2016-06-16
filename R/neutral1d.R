#' neutral1d
#'
#' Implements neutral model in one dimension
#'
#' @param size Size of the square grid on which to generate model. Total number
#' of points is size ^ 2
#' @param alpha Vector of two components respectively specifying the strength of
#' autocorrelation in time and space.
#' @param n Number of successive layers of temporal and spatial autocorrelation
#' used to generate final modelled values
#' @param sd0 Standard deviation of truncated normal distribution used to model
#' environmental variation (with mean of 1)
#' @param seed Random seed
#'
#' @return A vector of hotspot values sorted from high to low
#'
#' @seealso \code{neutral2d}
#'
#' @export
neutral1d <- function (size=10, alpha=c(0.1, 0.1), n=100, sd0=0.1, seed)
{
    if (alpha [1] <= 0)
        stop ('neutral model only makes sense with finite temporal autocorrelation')

    if (!missing (seed)) set.seed (seed)

    # generate truncated normal distribution in R and pass to Rcpp:
    yvec <- msm::rtnorm (size * size, mean=1, sd=sd0, lower=0, upper=2)

    y <- rcpp_neutral1d (size=size, alpha_t=alpha [1], alpha_s=alpha [2], 
                         nt=n, yvec=yvec)
    sort (y, decreasing=TRUE)
}
