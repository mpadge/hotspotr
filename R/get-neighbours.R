#' get_neighbours
#'
#' Determines neighbours of points (x, y) using Delaunay triangulation
#'
#' @param x Vector of x-coordinates
#' @param y Vector of y-coordinates
#'
#' @return A list of Delaunay neighbours of each point in (x, y)
#'
#' @examples
#' nbs <- get_neighbours (runif (10), runif (10))
#'
#' @export
get_neighbours <- function (x, y)
{
    if (!is.numeric (x) | !is.numeric (y))
        stop ('both x and y must be numeric')
    if (length (x) != length (y))
    {
        warning ('x and y should be same length; data will be trimmed')
        len <- min (c (length (x), length (y)))
        x <- x [seq (len)]
        y <- y [seq (len)]
    }

    rcpp_get_neighbours (x=x, y=y)
}
