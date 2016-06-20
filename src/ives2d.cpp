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

    Rcpp::NumericMatrix x (size, size), 
        xexp (size+2, size+2), xexp2 (size+2, size+2);
    std::fill (x.begin (), x.end (), 1.0);

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
        // spatial autocorrelation, calculated on expanded matrix
        for (int i=0; i<size; i++)
            for (int j=0; j<size; j++)
                xexp (i+1, j+1) = x (i, j);
        xexp (Rcpp::_, 0) = xexp (Rcpp::_, size);
        xexp (0, Rcpp::_) = xexp (size, Rcpp::_);
        xexp (Rcpp::_, size+1) = xexp (Rcpp::_, 1);
        xexp (size+1, Rcpp::_) = xexp (1, Rcpp::_);
        for (int i=0; i<size; i++)
            for (int j=0; j<size; j++)
                xexp2 (i+1, j+1) = (1.0 - 4.0 * alpha_s) * xexp (i+1, j+1) +
                    alpha_s * (xexp (i, j+1) + xexp (i+2, j+1) +
                            xexp (i+1, j) + xexp (i+1, j+2));
        x = xexp2 (Rcpp::Range (1, size), Rcpp::Range (1, size));
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

    Rcpp::NumericMatrix x (size, size), 
        xexp (size+2, size+2), xexp2 (size+2, size+2);
    std::fill (x.begin (), x.end (), 1.0);

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
        // spatial autocorrelation, calculated on expanded matrix
        for (int i=0; i<size; i++)
            for (int j=0; j<size; j++)
                xexp (i+1, j+1) = x (i, j);
        xexp (Rcpp::_, 0) = xexp (Rcpp::_, size);
        xexp (0, Rcpp::_) = xexp (size, Rcpp::_);
        xexp (Rcpp::_, size+1) = xexp (Rcpp::_, 1);
        xexp (size+1, Rcpp::_) = xexp (1, Rcpp::_);
        for (int i=0; i<size; i++)
            for (int j=0; j<size; j++)
            {
                vec = Rcpp::NumericVector::create (
                        xexp (i, j+1), xexp (i+2, j+1), 
                        xexp (i+1, j), xexp (i+1, j+2));
                tempd = Rcpp::max (vec);
                if (tempd < xexp (i+1, j+1))
                    tempd = Rcpp::min (vec);
                xexp2 (i+1, j+1) = (1.0 - alpha_s) * xexp (i+1, j+1) +
                    alpha_s * tempd;
            }
        x = xexp2 (Rcpp::Range (1, size), Rcpp::Range (1, size));
    }

    return x;
}
