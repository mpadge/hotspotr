#' morani
#'
#' Calculates Moran's I statistic for spatial association between 5 rectangular
#' neighbouring points on a regular grid
#'
#' @param z0 Square data matrix
#'
#' @return A vector of values of Moran's I statistic for each point in the grid
#'
#' @seealso \code{gearyc}, \code{getisord}
#'
#' @export
morani <- function (z0)
{
    size <- dim (z0) [1]
    z <- rbind (z0, z0, z0)
    z <- cbind (z, z, z)
    indx <- size + seq (size)
    
    xmn <- mean (z0)
    # local statistics
    wx <- (z [indx,indx] - xmn) * ((z [indx,indx] + z [indx-1,indx] +
                                    z[indx+1,indx] + z [indx,indx-1] +
                                    z[indx,indx+1]) - 4 * xmn) / 5
    s2 <- sum ((z0 - mean (z0)) ^ 2) / size ^ 2
    wx / s2 # global statistic is simply mean value
}
