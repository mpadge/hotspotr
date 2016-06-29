#include <Rcpp.h>
#include "ac-stats.h"
#include <fstream>

//' rcpp_neutral_hotspots
//'
//' Implements neutral model in two dimensions
//'
//' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
//' point. 
//' @param wts Weighting factors for each neighbour; must have same length as
//' nbs. Uniform weights used if not given.
//' @param alpha_t Strength of temporal autocorrelation
//' @param alpha_s Strength of spatial autocorrelation
//' @param sd0 Standard deviation of truncated normal distribution used to model
//' environmental variation (with mean of 1)
//' @param nt Number of successive layers of temporal and spatial autocorrelation
//' used to generate final modelled values
//'
//' @return A vector of simulated values of same size as \code{nbs}.
//'
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_neutral_hotspots (Rcpp::List nbs, Rcpp::List wts,
        double alpha_t, double alpha_s, double sd0, int nt)
{
    const int size = nbs.size ();

    int indx;
    double wtsum, tempd;

    // Set up truncated normal distribution
    Rcpp::NumericVector eps (nt * size);
    eps = Rcpp::rnorm (nt * size, 1.0, sd0);
    // TODO: Check timing with simple explicit loop instead of any
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

    Rcpp::NumericVector z (size), z2 (size);
    std::fill (z.begin (), z.end (), 1.0);

    Rcpp::NumericVector nbs1, wts1;
    for (int t=0; t<nt; t++)
    {
        // time step first
        for (int i=0; i<size; i++)
            z (i) = (1.0 - alpha_t) * z (i) + alpha_t * eps (t * size + i);
        // spatial autocorrelation
        z2 = Rcpp::clone (z);
        for (int i=0; i<size; i++)
        {
            tempd = wtsum = 0.0;
            nbs1 = Rcpp::as <Rcpp::NumericVector> (nbs (i));
            wts1 = Rcpp::as <Rcpp::NumericVector> (wts (i));
            for (int j=0; j<nbs1.size (); j++)
            {
                tempd += z (nbs1 (j) - 1) * wts1 (j);
                wtsum += wts1 (j);
            }
            z2 (i) = (1.0 - alpha_s) * z (i) + alpha_s * tempd / wtsum;
        }
        z = Rcpp::clone (z2);
    }

    return z;
}
