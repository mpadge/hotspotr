#include <Rcpp.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <vector>
/*
 * Examples of three different ways of extracting neighbourhood lists, and thus of
 * constructing adjacency matrices, from CGAL Delaunay triangulations. The first of
 * these is the easiest and neatest, but the others are very useful for illustrating
 * different ways of handing the CGAL routines.
 *
 * Has to be compiled with these linkages:
 * g++ aaa.cc -lCGAL -lgmp -frounding-math -o aaa
 */

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
int test ()
{
    // The second part of std::pair is extra info for the Delaunay points.
    std::vector< std::pair<Point,unsigned> > points;
    points.push_back( std::make_pair( Point(0.55,0.92), 100) );
    points.push_back( std::make_pair( Point(0.27,0.63), 101) );
    points.push_back( std::make_pair( Point(0.95,0.51), 102) );
    points.push_back( std::make_pair( Point(0.64,0.10), 103) );
    points.push_back( std::make_pair( Point(0.98,0.61), 104) );
    points.push_back( std::make_pair( Point(0.02,0.72), 105) );
    points.push_back( std::make_pair( Point(0.28,0.82), 106) );
    points.push_back( std::make_pair( Point(0.84,0.28), 107) );
    points.push_back( std::make_pair( Point(0.27,0.22), 108) );
    points.push_back( std::make_pair( Point(0.44,0.63), 109) );
          
    points.push_back( std::make_pair( Point(0.41,0.95), 110) );
    points.push_back( std::make_pair( Point(0.21,0.82), 111) );
    points.push_back( std::make_pair( Point(0.40,0.59), 112) );
    points.push_back( std::make_pair( Point(0.79,0.84), 113) );
    points.push_back( std::make_pair( Point(0.38,0.67), 114) );
    points.push_back( std::make_pair( Point(0.73,0.44), 115) );
    points.push_back( std::make_pair( Point(0.87,0.47), 116) );
    points.push_back( std::make_pair( Point(0.78,0.19), 117) );
    points.push_back( std::make_pair( Point(0.78,0.24), 118) );
    points.push_back( std::make_pair( Point(0.24,0.96), 119) );

    Delaunay triangulation;
    triangulation.insert (points.begin (), points.end ());

    /* This first way of extracting information uses a faces iterator to extract the
     * vertex->info for each face. A Delaunay triangulation has both finite and
     * infinite faces (see CGAL documentation). The infinite faces join to an
     * external, infinite vertex, so the finite_faces_iterator just includes the
     * internal faces. */
    for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin();
            fit != triangulation.finite_faces_end(); ++fit) {

        Delaunay::Face_handle face = fit;
        Rcpp::Rcout << "Triangle:\t" << triangulation.triangle(face) << std::endl;
        Rcpp::Rcout << "Vertex 0:\t" << triangulation.triangle(face)[0] << std::endl;
        Rcpp::Rcout << "[";
        for (int i=0; i<3; i++) 
        {
            Rcpp::Rcout << face->vertex(i)->info();
            if (i < 2) 
                Rcpp::Rcout << ", ";
            else 
                Rcpp::Rcout << "]" << std::endl;
        }
    }
    Rcpp::Rcout << "----------------------" << std::endl;

    // A vertex iterator can also be used:
    Vertex v;
    Vertex_iterator it = triangulation.vertices_begin (),
                        beyond = triangulation.vertices_end ();
    while (it != beyond) {
        v = *it;
        ++it;
        //Rcpp::Rcout << v.point () << std::endl;
    }
    //Rcpp::Rcout << "----------------------" << std::endl;

    /* And this is a long-hand way, using an iteration over all faces, checking
     * whether they are finite or not, and accessing the face.vertex info by
     * dereferencing pointers to each vertex of the face. */
    Face_iterator itf = triangulation.faces_begin (),
                  beyondf = triangulation.faces_end ();
    Face face;
    Face_handle neighbour;
    Fvertex_handle fvertex;

    Delaunay tr_reduced;
    int count = 0;
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
            Rcpp::Rcout << "(" << count << ") [ ";
            for (int i=0; i<3; i++) 
            {
                //fvertex = face.vertex (i);
                v = *(face.vertex (i));
                Rcpp::Rcout << v.info () << ", ";
            }
            Rcpp::Rcout << "]" << std::endl;
        }
        count++;
    }

    return 0;
}
