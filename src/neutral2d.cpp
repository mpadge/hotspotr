#include <Rcpp.h>
#include "ac_stats.h"

const double DOUBLE_MAX = std::numeric_limits<double>::max (),
    DOUBLE_MIN = -DOUBLE_MAX;

//' rcpp_neutral2d
//'
//' Implements neutral model in two dimensions
//'
//' @param size Size of the square grid on which to generate model. Total number
//' of points is size ^ 2
//' @param alpha_t Strength of temporal autocorrelation
//' @param alpha_s Strength of spatial autocorrelation
//' @param sd0 Standard deviation of truncated normal distribution used to model
//' environmental variation (with mean of 1)
//' @param nt Number of successive layers of temporal and spatial autocorrelation
//' used to generate final modelled values
//'
//' @return A matrix of dimension (size, size) of simulated values
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_neutral2d (int size, 
        double alpha_t, double alpha_s, double sd0, int nt)
{
    int indx;
    double tempd;

    // Set up truncated normal distribution
    Rcpp::NumericVector eps (nt * size * size);
    eps = Rcpp::rnorm (nt * size * size, 1.0, sd0);
    while (is_true (any (eps <= 0.0)))
    {
        tempd = Rcpp::rnorm (1, 1.0, sd0) (0);
        if (tempd >= 0.0 && tempd <= 2.0)
            eps (Rcpp::which_min (eps)) = tempd;
    }
    while (is_true (any (eps >= 2.0)))
    {
        tempd = Rcpp::rnorm (1, 1.0, sd0) (0);
        if (tempd >= 0.0 && tempd <= 2.0)
            eps (Rcpp::which_max (eps)) = tempd;
    }
        

    Rcpp::NumericMatrix x (size, size), x2 (size, size);
    std::fill (x.begin (), x.end (), 1.0);

    for (int t=0; t<nt; t++)
    {
        // temporal autocorrelation
        for (int i=1; i<size; i++)
            for (int j=0; j<size; j++)
            {
                indx = t * size * size + i * size + j;
                x (i, j) = (1.0 - alpha_t) * x (i, j) + alpha_t * eps (indx);
            }
        x2 = Rcpp::clone (x);
        // spatial autocorrelation
        for (int i=1; i<(size - 1); i++)
            for (int j=1; j<(size - 1); j++)
                x2 (i, j) = (1.0 - 4.0 * alpha_s) * x (i, j) + alpha_s *
                    (x (i-1, j) + x (i+1, j) + x (i, j-1) + x (i, j+1));
        // edges
        for (int i=1; i<(size - 1); i++)
        {
            x2 (i, 0) = (1.0 - 4.0 * alpha_s) * x (i, 0) + alpha_s *
                (x (i-1, 0) + x (i+1, 0) + x (i, size-1) + x (i, 1));
            x2 (i, size-1) = (1.0 - 4.0 * alpha_s) * x (i, size-1) + alpha_s *
                (x (i-1, size-1) + x (i+1, size-1) + x (i, size-2) + x (i, 0));
            x2 (0, i) = (1.0 - 4.0 * alpha_s) * x (0, i) + alpha_s *
                (x (size-1, i) + x (1, i) + x (0, i-1) + x (0, i+1));
            x2 (size-1, i) = (1.0 - 4.0 * alpha_s) * x (size-1, i) + alpha_s *
                (x (size-1, i-1) + x (size-1, i+1) + x (size-2, i) + x (0, i));
        }
        // corners
        x2 (0, 0) = (1.0 - 4.0 * alpha_s) * x (0, 0) + alpha_s *
            (x (0, 1) + x (0, size-1) + x (size-1, 0) + x (1, 0));
        x2 (0, size-1) = (1.0 - 4.0 * alpha_s) * x (0, size-1) + alpha_s *
            (x (0, size-2) + x (0, 0) + x (1, size-1) + x (size-1, size-1));
        x2 (size-1, 0) = (1.0 - 4.0 * alpha_s) * x (size-1, 0) + alpha_s *
            (x (size-1, 1) + x (size-1, size-1) + x (0, 0) + x (size-2, 0));
        x2 (size-1, size-1) = (1.0 - 4.0 * alpha_s) * x (size-1, size-1) + 
                alpha_s * (x (size-2, size-1) + x (0, size-1) + 
                x (size-1, size-2) + x (size-1, 0));

        x = Rcpp::clone (x2);
    } // end for t

    return x;
}


//' rcpp_neutral2d_ntests
//'
//' Performs repeated neutral tests to yield average distributions of both
//' hotspot values and spatial autocorrelation statistics.
//'
//' @param size Size of the square grid on which to generate model. Total number
//' of points is size ^ 2
//' @param alpha_t Strength of temporal autocorrelation
//' @param alpha_s Strength of spatial autocorrelation
//' @param sd0 Standard deviation of truncated normal distribution used to model
//' environmental variation (with mean of 1)
//' @param nt Number of successive layers of temporal and spatial autocorrelation
//' used to generate final modelled values
//' @param ntests Number of tests used to obtain average values
//' @param ac_type Character string specifying type of aucorrelation
//' (\code{moran}, \code{geary}, or code{getis-ord}).
//'
//' @return A matrix of dimension (size * size, 2), with first column containing
//' sorted and re-scaled hotspot values, and second column containing sorted and
//' re-scaled spatial autocorrelation statistics.
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_neutral2d_ntests (int size, 
        double alpha_t, double alpha_s, double sd0, int nt, int ntests,
        std::string ac_type)
{
    // TODO Change int ac_type to std::string
    int len = size * size;
    double xmin, xmax;

    Rcpp::NumericMatrix x (size, size);

    Rcpp::NumericVector xv (len), xv_tot (len), ac_vec (len), ac_vec1 (len);
    std::fill (ac_vec.begin (), ac_vec.end (), 0.0);
    std::fill (xv_tot.begin (), xv_tot.end (), 0.0);
    for (int n=0; n<ntests; n++)
    {
        x = rcpp_neutral2d (size, alpha_t, alpha_s, sd0, nt);
        ac_vec1 = rcpp_ac_stats (x, ac_type); // already sorted and normalised
        for (int i=0; i<size; i++)
            for (int j=0; j<size; j++)
                xv (i * size + j) = x (i, j);
        for (int i=0; i<len; i++)
            ac_vec (i) = ac_vec (i) + ac_vec1 (i);
        std::sort (xv.begin (), xv.end (), std::greater<double> ());
        xmin = Rcpp::min (xv);
        xmax = Rcpp::max (xv);
        for (int i=0; i<len; i++)
            xv_tot (i) = xv_tot (i) + (xv (i) - xmin) / (xmax - xmin);
    }

    Rcpp::NumericMatrix result (len, 2);
    for (int i=0; i<len; i++)
    {
        result (i, 0) = xv_tot (i) / (double) ntests;
        result (i, 1) = ac_vec (i) / (double) ntests;
    }
    Rcpp::colnames (result) = Rcpp::CharacterVector::create ("y", "ac");

    return result;
}
