#' gearyc
#'
#' Calculates Geary's C statistic for spatial association between 5 rectangular
#' neighbouring points on a regular grid
#'
#' @param z0 Square data matrix
#'
#' @return A vector of values of Geary's C statistic for each point in the grid
#'
#' @seealso \code{morani}, \code{getisord}
#'
#' @export
gearyc <- function (z0)
{
    size <- dim (z0) [1]
    z <- rbind (z0, z0, z0)
    z <- cbind (z, z, z)
    indx <- size + seq (size)
    
    xmn <- mean (z0)
    wx <- (z [indx,indx] - z [indx-1,indx]) ^ 2 +
            (z [indx,indx] - z [indx+1,indx]) ^ 2 +
            (z [indx,indx] - z [indx,indx-1]) ^ 2 +
            (z [indx,indx] - z [indx,indx+1]) ^ 2
    wx <- wx / 5
    s2 <- sum ((z0 - mean (z0)) ^ 2) / size ^ 2
    wx / (2 * s2) # global statistic is simply mean value
}
