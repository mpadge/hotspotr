#' neutral2d
#'
#' Implements neutral model in two dimensions
#'
#' @param size Size of the square grid on which to generate model. Total number
#' of points is size ^ 2
#' @param alpha Vector of two components respectively specifying the strength of
#' autocorrelation in time and space.
#' @param n Number of successive layers of temporal and spatial autocorrelation
#' used to generate final modelled values
#' @param seed Random seed
#'
#' @return A vector of hotspot values sorted from high to low
#'
#' @seealso \code{neutral1d}
#'
#' @export
neutral2d <- function (size=10, alpha=c(0.1, 0.1), n=100, seed)
{
    if (!missing (seed)) set.seed (seed)

    yn <- rep (1, size * size)
    for (i in 1:n)
    {
        # temporal autocorrelation with alpha [1]
        yn <- yn * rnorm (size * size, 1, 0.1)
        indx <- 1:(length (yn) - 1)
        yn [indx+1] <- (1 - alpha [1]) * yn [indx] + alpha [1] * yn [indx + 1]
        # spatial autocorrelation with alpha [2]
        yn <- matrix (yn, nrow=size, ncol=size)
        yn2 <- cbind (yn [,size], yn, yn [,1])
        yn2 <- rbind (yn2 [size,], yn2, yn2 [1,])
        indx <- 1:size + 1
        yn2 [indx, indx] <- (1 - 4 * alpha [2]) * yn2 [indx, indx] + alpha [2] *
                    (yn2 [indx-1, indx] + yn2 [indx+1, indx] +
                     yn2 [indx, indx-1] + yn2 [indx, indx+1])
        yn <- yn2 [indx, indx]
    }

    sort (yn, decreasing=TRUE)
}
