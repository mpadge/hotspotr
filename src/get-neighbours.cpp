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

typedef Tds::Vertex Vertex;
typedef Tds::Vertex_iterator Vertex_iterator;
typedef Tds::Face Face;
typedef Tds::Face_iterator Face_iterator;
typedef Face::Face_handle Face_handle;
typedef Face::Vertex_handle Fvertex_handle;

// [[Rcpp::export]]
Rcpp::List rcpp_get_neighbours (Rcpp::NumericVector x, Rcpp::NumericVector y)
{
    // The second part of std::pair is extra info for the Delaunay points.
    std::vector< std::pair<Point, unsigned> > points;
    unsigned number = 0;
    for (int i=0; i<x.size (); i++)
    {
        points.push_back ( std::make_pair (Point (x (i), y (i)), ++number) );
    }

    Delaunay triangulation;
    triangulation.insert (points.begin (), points.end ());

    /* This first way of extracting information uses a faces iterator to extract the
     * vertex->info for each face. A Delaunay triangulation has both finite and
     * infinite faces (see CGAL documentation). The infinite faces join to an
     * external, infinite vertex, so the finite_faces_iterator just includes the
     * internal faces. */
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

    /* And this is a long-hand way, using an iteration over all faces, checking
     * whether they are finite or not, and accessing the face.vertex info by
     * dereferencing pointers to each vertex of the face. */
    Vertex v;
    Face_iterator itf = triangulation.faces_begin (),
                  beyondf = triangulation.faces_end ();
    Face face;
    Face_handle neighbour;
    Fvertex_handle fvertex;

    Delaunay tr_reduced;
    bool is_finite;
    while (itf != beyondf) 
    {
        face = *(itf++);
        is_finite = true;
        for (int i=0; i<3; i++) {
            fvertex = face.vertex (i);
            if (triangulation.is_infinite (fvertex)) 
                is_finite = false;
        }
        if (is_finite)
        {
            for (int i=0; i<3; i++) 
                v = *(face.vertex (i));
        }
    }

    points.resize (0);

    return nbs;
}
