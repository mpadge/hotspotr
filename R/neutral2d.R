#' neutral2d
#'
#' Implements neutral model in two dimensions
#'
#' @param size Size of the square grid on which to generate model. Total number
#' of points is size ^ 2
#' @param alpha Vector of two components respectively specifying the strength of
#' autocorrelation in time and space.
#' @param nt Number of successive layers of temporal and spatial autocorrelation
#' used to generate final modelled values
#' @param sd0 Standard deviation of truncated normal distribution used to model
#' environmental variation (with mean of 1)
#' @param seed Random seed
#'
#' @return A vector of hotspot values sorted from high to low
#'
#' @seealso \code{neutral1d}
#'
#' @examples
#' y <- neutral2d ()
#'
#' @export
neutral2d <- function (size=10, alpha=c(0.1, 0.1), nt=100, sd0=0.1, seed)
{
    if (alpha [1] <= 0)
        stop ('neutral model only makes sense with finite temporal autocorrelation')

    if (!missing (seed)) set.seed (seed)

    eps <- rnorm (size * size * nt, mean=0, sd=sd0)
    y <- rcpp_neutral2d (size=size, alpha_t=alpha [1], alpha_s=alpha [2], 
                         nt=nt, eps=eps)
    matrix (y, nrow=size, ncol=size)
}
