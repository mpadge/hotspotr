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
#'
#' @return A vector of hotspot values sorted from high to low
#'
#' @seealso \code{neutral2d}
#'
#' @export
neutral1d <- function (size=10, alpha=c(0.1, 0.1), n=1)
{
    yn <- rep (1, size * size)
    for (i in 1:n)
    {
        # temporal autocorrelation with alpha [1]
        yn <- yn * rnorm (size * size, 1, 0.1)
        indx <- 1:(length (yn) - 1)
        yn [indx+1] <- (1 - alpha [1]) * yn [indx] + alpha [1] * yn [indx + 1]
        # spatial autocorrelation with alpha [2]
        yn <- c (yn [size * size], yn, yn [1])
        indx <- 1:(size^2) + 1
        yn [indx] <- (1 - 2 * alpha [2]) * yn [indx] +
                        alpha [2] * (yn [indx - 1] + yn [indx + 1])
        yn <- yn [indx]
    }

    sort (yn, decreasing=TRUE)
}
