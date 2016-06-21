#include <Rcpp.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel            Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                       Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>                    Delaunay;
typedef Kernel::Point_2                                                Point;

//' rcpp_get_neighbours
//'
//' Determines neighbours of points (x, y) using Delaunay triangulation
//'
//' @param x Vector of x-coordinates
//' @param y Vector of y-coordinates
//'
//' @return List of Delaunay triangle membership
//'
// [[Rcpp::export]]
Rcpp::List rcpp_get_neighbours (Rcpp::NumericVector x, Rcpp::NumericVector y)
{
    // The second part of std::pair is extra info for the Delaunay points.
    std::vector< std::pair<Point, unsigned> > points;
    unsigned number = 0;
    for (int i=0; i<x.size (); i++)
        points.push_back ( std::make_pair (Point (x (i), y (i)), number++) );

    Delaunay triangulation;
    triangulation.insert (points.begin (), points.end ());

    Rcpp::List nbs;

    int ntri = 0;
    for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin();
            fit != triangulation.finite_faces_end(); ++fit) 
        ntri++;

    Rcpp::NumericMatrix triangles (ntri, 3);
    ntri = 0;
    for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin();
            fit != triangulation.finite_faces_end(); ++fit) 
    {
        Delaunay::Face_handle face = fit;
        for (int i=0; i<3; i++) 
            triangles (ntri, i) = face -> vertex (i) -> info ();
        ntri++;
    }

    /* These triangles are then converted to lists of neighbouring points. The
     * following code is easily implemented in R, assuming rcpp_get_neighbours
     * just returns a list of triangular neighbours (that is, triangles in the
     * above constructed as an Rcpp::List instead of Rcpp::Matrix):
     * > nbs <- do.call (rbind, rcpp_get_neighbours (x=x, y=y))
     * > nbs2 <- list ()
     * > for (i in seq (length (x)))
     * > {
     * >     indx <- which (apply (nbs, 1, function (j) i %in%
     * > )
     * >     pts <- sort (unique (as.numeric (nbs [indx,])))
     * >     nbs2 [[i]] <- pts [!pts %in% i]
     * > }
     * ----------------- microbenchmark results:
     *  R = 62.8ms, Rcpp = 0.34ms
     * The Rcpp speed-up is over 99.5%!!!
     */

    int tempi;
    std::vector <int> nbsi;
    bool i_in_list;

    for (int i=0; i<x.size (); i++)
    {
        nbsi.resize (0);
        for (int j=0; j<triangles.nrow (); j++)
        {
            i_in_list = false;
            for (int k=0; k<3; k++)
                if (triangles (j, k) == i)
                    i_in_list = true;
            if (i_in_list)
                for (int k=0; k<3; k++)
                {
                    tempi = triangles (j, k);
                    if (tempi != i && std::find (nbsi.begin (), nbsi.end (),
                                tempi) == nbsi.end ())
                        nbsi.push_back (tempi);
                }
        } // end for j over nrow (triangles)
        std::sort (nbsi.begin (), nbsi.end ());
        // Then adjust nbsi by 1 to make them R-indexed
        for (std::vector <int>::iterator it=nbsi.begin (); 
                it != nbsi.end (); ++it)
            *it = *it + 1;
        nbs.push_back (nbsi);
    }

    nbsi.resize (0);
    points.resize (0);

    return nbs;
}
