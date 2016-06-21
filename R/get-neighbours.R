#' get_neighbours
#'
#' Determines neighbours of points (x, y) using Delaunay triangulation
#'
#' @param x Vector of x-coordinates
#' @param y Vector of y-coordinates
#' @param plot Produce a plot of the Delaunay triangulation
#'
#' @return A vector of hotspot values sorted from high to low
#'
#' @examples
#' nbs <- get_neighbours (runif (10), runif (10))
#'
#' @export
get_neighbours <- function (x, y, plot=TRUE)
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

    nbs <- do.call (rbind, rcpp_get_neighbours (x=x, y=y))
    if (plot)
    {
        plot (x, y, pch=1, xlab="", ylab="")
        cols <- rainbow (nrow (nbs))
        ps0 <- par ('ps')
        par (ps=20, font=2)
        for (i in seq (nrow (nbs)))
        {
            indx <- c (nbs [i,], nbs [i,1])
            lines (x [indx], y [indx], col=cols [i])
            text (mean (x [indx [1:3]]), mean (y [indx [1:3]]),
                  col=cols [i], labels=i)
        }
        par (ps=16, font=1)
        text (x, y, labels=seq (length (x)))
        par (ps=ps0)
    }
    # nbs is a list of triangular neighbours for each point. This is converted
    # to list of *all* neighbours for each point
    nbs2 <- list ()
    for (i in seq (length (x)))
    {
        indx <- which (apply (nbs, 1, function (j) i %in% j))
        pts <- sort (unique (as.numeric (nbs [indx,])))
        nbs2 [[i]] <- pts [!pts %in% i]
    }
    return (nbs2)
}
