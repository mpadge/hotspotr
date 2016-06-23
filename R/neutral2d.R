#' neutral2d
#'
#' Implements neutral model in two dimensions
#'
#' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
#' point
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
neutral2d <- function (nbs, alpha=c(0.1, 0.1), nt=100, sd0=0.1, seed)
{
    if (alpha [1] <= 0)
        stop ('neutral model only makes sense with finite temporal autocorrelation')

    if (!missing (seed)) set.seed (seed)

    rcpp_neutral2d (nbs=nbs, alpha_t=alpha [1], alpha_s=alpha [2],
                    sd0=sd0, nt=nt)
}
