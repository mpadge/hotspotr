#' ives2d
#'
#' Generates a 2D surface using the temporal autocorrelation model of Ives &
#' Klopfer (Ecology 1997).
#'
#' @param size Size of the square grid on which to generate model. Total number
#' of points is size ^ 2
#' @param nt Number of successive layers of temporal and spatial autocorrelation
#' used to generate final modelled values
#' @param sd0 Standard deviation of distributions from which values of 'r' and
#' 's' are drawn.
#' @param alpha Vector of two components respectively specifying the strength of
#' autocorrelation in time and space.
#' @param spatial If TRUE, spatial autocorrelation is calculated by assuming
#' movement along locally maximal gradients, otherwise it is equivalent to
#' neutral model
#' @param seed Random seed
#'
#' @return A matrix of (size, size)
#'
#' @export
ives2d <- function (size=10, nt=1000, sd0=0.1, alpha=c(0.1, 0.1), spatial=FALSE,
                    seed) 
{
    if (!missing (seed)) set.seed (seed)
    s0 <- 0.5
    r0 <- 1.05
    # generate truncated normal distributions in R and pass them to Rcpp:
    svec <- msm::rtnorm (nt * size * size, mean=s0, sd=sd0, lower=0, upper=2*s0)
    rvec <- msm::rtnorm (nt * size * size, mean=r0, sd=sd0, lower=0, upper=2*r0)

    result <- NULL
    if (!spatial)
        result <- rcpp_ives2d (size, nt, alpha[1], alpha[2], svec, rvec)
    else
        result <- rcpp_ives2d_space (size, nt, alpha[1], alpha[2], svec, rvec)
    return (result)
}
