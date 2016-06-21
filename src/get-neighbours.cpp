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

// [[Rcpp::export]]
Rcpp::List rcpp_get_neighbours (Rcpp::NumericVector x, Rcpp::NumericVector y)
{
    // The second part of std::pair is extra info for the Delaunay points.
    std::vector< std::pair<Point, unsigned> > points;
    unsigned number = 0;
    for (int i=0; i<x.size (); i++)
        points.push_back ( std::make_pair (Point (x (i), y (i)), ++number) );

    Delaunay triangulation;
    triangulation.insert (points.begin (), points.end ());

    Rcpp::List nbs;
    std::vector <int> nbsi;
    nbsi.resize (0);
    for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin();
            fit != triangulation.finite_faces_end(); ++fit) 
    {
        nbsi.resize (0);
        Delaunay::Face_handle face = fit;
        for (int i=0; i<3; i++) 
            nbsi.push_back (face -> vertex (i) -> info ());
        nbs.push_back (nbsi);
    }
    nbsi.resize (0);

    points.resize (0);

    return nbs;
}
