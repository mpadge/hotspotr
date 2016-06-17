#' morani
#'
#' Calculates Moran's I statistic for spatial association between 5 rectangular
#' neighbouring points on a regular grid
#'
#' @param z0 Square data matrix
#'
#' @return A sorted vector of values of Moran's I statistic re-scaled between 0
#' and 1
#'
#' @seealso \code{gearyc}, \code{getisord}
#'
#' @export
morani <- function (z0)
{
    size <- dim (z0) [1]

    z <- rbind (z0 [size,], z0, z0 [1,])
    z <- cbind (z [,size], z, z [,1])
    indx <- seq (size) + 1
    
    xmn <- mean (z0)
    # local statistics
    wx <- (z [indx,indx] - xmn) * ((z [indx,indx] + z [indx-1,indx] +
                                    z[indx+1,indx] + z [indx,indx-1] +
                                    z[indx,indx+1]) - 4 * xmn) / 5
    s2 <- sum ((z0 - mean (z0)) ^ 2) / size ^ 2
    wx <- wx / s2 # global statistic is simply mean value
    sort ((wx - min (wx)) / diff (range (wx)), decreasing=TRUE)
}
