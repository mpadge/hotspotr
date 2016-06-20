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
        

    Rcpp::NumericMatrix x (size, size), 
        xexp (size+2, size+2), xexp2 (size+2, size+2);
    std::fill (x.begin (), x.end (), 1.0);

    for (int t=0; t<nt; t++)
    {
        // temporal autocorrelation
        for (int i=0; i<size; i++)
            for (int j=0; j<size; j++)
            {
                indx = t * size * size + i * size + j;
                x (i, j) = (1.0 - alpha_t) * x (i, j) + alpha_t * eps (indx);
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
