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
            //z2 (i) = (1.0 - (double) nbs1.size () * alpha_s) * z (i) +
            //    alpha_s * (double) nbs1.size() * tempd / wtsum;
            z2 (i) = (1.0 - alpha_s) * z (i) + alpha_s * tempd / wtsum;
        }
        z = Rcpp::clone (z2);
    }

    return z;
}


//' rcpp_neutral_hotspots_ntests
//'
//' Performs repeated neutral tests to yield average distributions of both
//' hotspot values and spatial autocorrelation statistics.
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
//' @param ntests Number of tests used to obtain average values
//' @param ac_type Character string specifying type of aucorrelation
//' (\code{moran}, \code{geary}, or code{getis-ord}).
//'
//' @return A matrix of dimension (size, 2), with first column containing
//' sorted and re-scaled hotspot values, and second column containing sorted and
//' re-scaled spatial autocorrelation statistics.
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_neutral_hotspots_ntests (Rcpp::List nbs, 
        Rcpp::List wts, double alpha_t, double alpha_s, double sd0, int nt, int
        ntests, std::string ac_type)
{
    const int size = nbs.size ();

    Rcpp::NumericVector z (size), z1 (size), ac (size), ac1 (size);
    std::fill (ac.begin (), ac.end (), 0.0);
    std::fill (z.begin (), z.end (), 0.0);

    for (int n=0; n<ntests; n++)
    {
        z1 = rcpp_neutral_hotspots (nbs, wts, alpha_t, alpha_s, sd0, nt);
        ac1 = rcpp_ac_stats (z1, nbs, wts, ac_type); // sorted and normalised
        ac += ac1;
        std::sort (z1.begin (), z1.end (), std::greater<double> ());
        z += (z1 - (double) Rcpp::min (z1)) /
                ((double) Rcpp::max (z1) - (double) Rcpp::min (z1));
    }

    Rcpp::NumericMatrix result (size, 2);
    result (Rcpp::_, 0) = z / (double) ntests;
    result (Rcpp::_, 1) = ac / (double) ntests;
    Rcpp::colnames (result) = Rcpp::CharacterVector::create ("y", "ac");

    return result;
}
