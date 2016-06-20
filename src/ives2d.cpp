#include <Rcpp.h>
#include <fstream>
#include <algorithm> // std::max

//' rcpp_ives2d
//'
//' Implements neutral model in two dimensions
//'
//' @param size Size of the square grid on which to generate model. Total number
//' of points is size ^ 2
//' @param nt Number of successive layers of temporal and spatial autocorrelation
//' used to generate final modelled values
//' @param alpha_t Strength of temporal autocorrelation
//' @param alpha_s Strength of spatial autocorrelation
//' @param svec Vector of random numbers for values of \code{s} drawn from a
//' truncated normal ' distribution.
//' @param rvec Vector of random numbers for values of \code{r} drawn from a
//' truncated normal ' distribution.
//'
//' @return A matrix of dimension (size, size) of simulated values
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_ives2d (int size, int nt, 
        double alpha_t, double alpha_s,
        Rcpp::NumericVector svec, Rcpp::NumericVector rvec)
{
    int indx;
    const double a0 = 0.05;

    if (svec.size () != (size * size * nt))
        Rcpp::stop ("svec must have length (size * size * nt)");
    if (rvec.size () != (size * size * nt))
        Rcpp::stop ("rvec must have length (size * size * nt)");

    Rcpp::NumericMatrix x (size, size), x2 (size, size);
    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
            x (i, j) = 1.0;

    for (int t=0; t<nt; t++)
    {
        // time step first
        for (int i=0; i<size; i++)
            for (int j=0; j<size; j++)
            {
                indx = t * size * size + i * size + j;
                x (i, j) = svec [indx] * x (i, j) * 
                    (1.0 + rvec [indx] / (1.0 + a0 * x(i, j)));
            }
        x2 = Rcpp::clone (x);
        // then spatial autocorrelation
        for (int i=1; i<(size - 1); i++)
            for (int j=1; j<(size - 1); j++)
                x2 (i, j) = (1.0 - 4.0 * alpha_s) * x (i ,j) + alpha_s * 
                    (x(i-1, j) + x(i+1, j) + x(i, j-1) + x(i, j+1));
        // wrap edges
        for (int i=1; i<(size - 1); i++)
        {
            x2 (i, 0) = (1.0 - 4.0 * alpha_s) * x (i ,0) + alpha_s * 
                (x(i-1, 0) + x(i+1, 0) + x(i, size-1) + x(i, 1));
            x2 (i, size-1) = (1.0 - 4.0 * alpha_s) * x (i ,size-1) + alpha_s * 
                (x(i-1, size-1) + x(i+1, size-1) + x(i, size-2) + x(i, 0));
            x2 (0, i) = (1.0 - 4.0 * alpha_s) * x (0 ,i) + alpha_s * 
                (x(0, i-1) + x(0, i+1) + x(size-1, i) + x(1, i));
            x2 (size-1, i) = (1.0 - 4.0 * alpha_s) * x (size-1, i) + alpha_s * 
                (x(size-1, i-1) + x(size-1, i+1) + x(size-2, i) + x(0, i));
        }
        // and corners
        x2 (0, 0) = (1.0 - 4.0 * alpha_s) * x (0, 0) + alpha_s *
            (x (0, 1) + x (0, size-1) + x (size-1, 0) + x (1, 0));
        x2 (0, size-1) = (1.0 - 4.0 * alpha_s) * x (0, size-1) + alpha_s *
            (x (0, size-2) + x (0, 0) + x (1, size-1) + x (size-1, size-1));
        x2 (size-1, 0) = (1.0 - 4.0 * alpha_s) * x (size-1, 0) + alpha_s *
            (x (size-1, 1) + x (size-1, size-1) + x (0, 0) + x (size-2, 0));
        x2 (size-1, size-1) = (1.0 - 4.0 * alpha_s) * x (size-1, size-1) + 
            alpha_s * (x (size-2, size-1) + x (0, size-1) + 
                    x (size-1, 0) + x (size-1, size-2));

        x = Rcpp::clone (x2);
    }

    return x;
}
//' rcpp_ives2d_space
//'
//' Implements neutral model in two dimensions with additional spatial
//' structure, implented here through replacing generic local autocorrelation
//' with movement along maximal local gradients.
//'
//' @param size Size of the square grid on which to generate model. Total number
//' of points is size ^ 2
//' @param nt Number of successive layers of temporal and spatial autocorrelation
//' used to generate final modelled values
//' @param alpha_t Strength of temporal autocorrelation
//' @param alpha_s Strength of spatial autocorrelation
//' @param svec Vector of random numbers for values of \code{s} drawn from a
//' truncated normal ' distribution.
//' @param rvec Vector of random numbers for values of \code{r} drawn from a
//' truncated normal ' distribution.
//'
//' @return A matrix of dimension (size, size) of simulated values
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_ives2d_space (int size, int nt, 
        double alpha_t, double alpha_s,
        Rcpp::NumericVector svec, Rcpp::NumericVector rvec)
{
    const double a0 = 0.05;

    int indx;
    double tempd;

    if (svec.size () != (size * size * nt))
        Rcpp::stop ("svec must have length (size * size * nt)");
    if (rvec.size () != (size * size * nt))
        Rcpp::stop ("rvec must have length (size * size * nt)");

    Rcpp::NumericMatrix x (size, size), x2 (size, size);
    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
            x (i, j) = 1.0;

    Rcpp::NumericVector vec;
    for (int t=0; t<nt; t++)
    {
        // time step first
        for (int i=0; i<size; i++)
            for (int j=0; j<size; j++)
            {
                indx = t * size * size + i * size + j;
                x (i, j) = svec [indx] * x (i, j) * 
                    (1.0 + rvec [indx] / (1.0 + a0 * x(i, j)));
            }
        x2 = Rcpp::clone (x);
        // then spatial autocorrelation
        for (int i=1; i<(size - 1); i++)
            for (int j=1; j<(size - 1); j++)
            {
                vec = Rcpp::NumericVector::create (
                        x(i-1, j), x(i+1, j), x(i, j-1), x(i, j+1));
                tempd = Rcpp::max (vec);
                if (tempd < x (i, j))
                    tempd = Rcpp::min (vec);
                x2 (i, j) = (1.0 - alpha_s) * x (i, j) + alpha_s * tempd;
            }
        // wrap edges
        for (int i=1; i<(size - 1); i++)
        {
            vec = Rcpp::NumericVector::create (
                    x(i-1, 0), x(i+1, 0), x(i, size-1), x(i, 1));
            tempd = Rcpp::max (vec);
            if (tempd < x (i, 0))
                tempd = Rcpp::min (vec);
            x2 (i, 0) = (1.0 - alpha_s) * x (i, 0) + alpha_s * tempd;

            vec = Rcpp::NumericVector::create (
                    x(i-1, size-1), x(i+1, size-1), x(i, size-2), x(i, 0));
            tempd = Rcpp::max (vec);
            if (tempd < x (i, size-1))
                tempd = Rcpp::min (vec);
            x2 (i, size-1) = (1.0 - alpha_s) * x (i, size-1) + alpha_s * tempd;

            vec = Rcpp::NumericVector::create (
                    x(0, i-1), x(0, i+1), x(size-1, i), x(1, i));
            tempd = Rcpp::max (vec);
            if (tempd < x (0, i))
                tempd = Rcpp::min (vec);
            x2 (0, i) = (1.0 - alpha_s) * x (0, i) + alpha_s * tempd;

            vec = Rcpp::NumericVector::create (
                    x(size-1, i-1), x(size-1, i+1), x(size-2, i), x(0, i));
            tempd = Rcpp::max (vec);
            if (tempd < x (size-1, i))
                tempd = Rcpp::min (vec);
            x2 (size-1, i) = (1.0 - alpha_s) * x (size-1, i) + alpha_s * tempd;
        }
        // and corners
        vec = Rcpp::NumericVector::create (
                x (0, 1), x (0, size-1), x (size-1, 0), x (1, 0));
        tempd = Rcpp::max (vec);
        if (tempd < x (0, 0))
            tempd = Rcpp::min (vec);
        x2 (0, 0) = (1.0 - alpha_s) * x (0, 0) + alpha_s * tempd;

        vec = Rcpp::NumericVector::create (
                x (0, 1), x (0, size-1), x (size-1, 0), x (1, 0));
        tempd = Rcpp::max (vec);
        if (tempd < x (0, size-1))
            tempd = Rcpp::min (vec);
        x2 (0, size-1) = (1.0 - alpha_s) * x (0, size-1) + alpha_s * tempd;

        x2 (0, size-1) = (1.0 - 4.0 * alpha_s) * x (0, size-1) + alpha_s *
            (x (0, size-2) + x (0, 0) + x (1, size-1) + x (size-1, size-1));
        x2 (size-1, 0) = (1.0 - 4.0 * alpha_s) * x (size-1, 0) + alpha_s *
            (x (size-1, 1) + x (size-1, size-1) + x (0, 0) + x (size-2, 0));
        x2 (size-1, size-1) = (1.0 - 4.0 * alpha_s) * x (size-1, size-1) + 
            alpha_s * (x (size-2, size-1) + x (0, size-1) + 
                    x (size-1, 0) + x (size-1, size-2));

        x = Rcpp::clone (x2);
    }

    return x;
}
