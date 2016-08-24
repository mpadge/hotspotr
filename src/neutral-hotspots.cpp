#include <Rcpp.h>
#include "ac-stats.h"
#include <fstream>

//' rcpp_trunc_ndist
//'
//' Truncated normal distribution (mean 1, respective upper and lower limits of
//' 0 and 2).
//'
//' @param len Number of elements to be simulated
//' @param sd Standard deviation
//'
//' @return A vector of truncated normally distributed values
//'
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_trunc_ndist (int len, double sd)
{
    double tempd;

    // Set up truncated normal distribution
    Rcpp::NumericVector eps (len);
    eps = Rcpp::rnorm (len, 1.0, sd);
    // TODO: Check timing with simple explicit loop instead of any
    while (Rcpp::is_true (Rcpp::any (eps <= 0.0)))
    {
        tempd = Rcpp::rnorm (1, 1.0, sd) (0);
        if (tempd >= 0.0 && tempd <= 2.0)
            eps (Rcpp::which_min (eps)) = tempd;
    }
    while (Rcpp::is_true (Rcpp::any (eps >= 2.0)))
    {
        tempd = Rcpp::rnorm (1, 1.0, sd) (0);
        if (tempd >= 0.0 && tempd <= 2.0)
            eps (Rcpp::which_min (eps)) = tempd;
    }

    return eps;
}

//' rcpp_neutral_hotspots
//'
//' Implements neutral model in two dimensions
//'
//' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
//' point. 
//' @param wts Weighting factors for each neighbour; must have same length as
//' nbs. 
//' @param nbsi List of matrices as returned from \code{get_nbsi}. each element
//' of which contains the i-th nearest neighbour to each point.
//' @param alpha Strength of spatial autocorrelation
//' @param sd0 Standard deviation of truncated normal distribution used to model
//' environmental variation (with mean of 1)
//' @param log_scale If TRUE, raw hotspot values are log-transformed
//' @param niters Number of iterations of spatial autocorrelation
//' @param ac_type Character string specifying type of aucorrelation
//' (\code{moran}, \code{geary}, or code{getis-ord}).
//'
//' @return A vector of simulated values of same size as \code{nbs}.
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_neutral_hotspots (Rcpp::List nbs, Rcpp::List wts,
        Rcpp::List nbsi, double alpha, double sd0, bool log_scale, int niters,
        std::string ac_type)
{
    const int size = nbs.size ();

    int indx;
    double wtsum, tempd;

    Rcpp::NumericMatrix tempmat;
    Rcpp::NumericVector z = rcpp_trunc_ndist (size, sd0); 
    Rcpp::NumericVector z_copy (size), nbs_to, nbs_from, nbs_n, ac;

    // Spatial autocorrelation
    for (int i=0; i<niters; i++)
    {
        std::fill (z_copy.begin (), z_copy.end (), 0.0);
        for (int j=0; j<nbsi.size (); j++)
        {
            tempmat = Rcpp::as <Rcpp::NumericMatrix> (nbsi (j));
            nbs_to = tempmat (Rcpp::_, 0);
            nbs_from = tempmat (Rcpp::_, 1);
            nbs_n = tempmat (Rcpp::_, 2);
            // Note that nbs are 1-indexed
            for (int k=0; k<nbs_to.size (); k++)
            {
                z_copy (nbs_to (k) - 1) += (1.0 - alpha) * z (nbs_to (k) - 1) +
                    alpha * z (nbs_from (k) - 1) / (double) nbs_n (k);
            }
            std::copy (z.begin (), z.end (), z_copy.begin ());
        } // end for j over nbsi
    } // end for i over niters

    if (log_scale)
        z = log10 (z);

    ac = rcpp_ac_stats (z, nbs, wts, ac_type);
    std::sort (z.begin (), z.end (), std::greater<double> ());
    
    z = (z - (double) Rcpp::min (z)) /
        ((double) Rcpp::max (z) - (double) Rcpp::min (z));

    Rcpp::NumericMatrix result (size, 2);
    result (Rcpp::_, 0) = z;
    result (Rcpp::_, 1) = ac;
    Rcpp::colnames (result) = Rcpp::CharacterVector::create ("y", "ac");

    return result;
}

//' rcpp_neutral_hotspots_ntests
//'
//' Performs repeated neutral tests to yield average distributions of both
//' hotspot values and spatial autocorrelation statistics.
//'
//' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
//' point. 
//' @param wts Weighting factors for each neighbour; must have same length as
//' nbs. 
//' @param nbsi List of matrices as returned from \code{get_nbsi}. each element
//' of which contains the i-th nearest neighbour to each point.
//' @param alpha Strength of spatial autocorrelation
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
        Rcpp::List wts, Rcpp::List nbsi, double alpha, double sd0, int niters,
        std::string ac_type, bool log_scale, int ntests)
{
    const int size = nbs.size ();

    Rcpp::NumericMatrix hs1;
    Rcpp::NumericVector z (size), z1 (size), ac (size), ac1 (size);
    std::fill (ac.begin (), ac.end (), 0.0);
    std::fill (z.begin (), z.end (), 0.0);

    for (int n=0; n<ntests; n++)
    {
        hs1 = rcpp_neutral_hotspots (nbs, wts, nbsi, alpha, sd0, log_scale,
                niters, ac_type);
        z += hs1 (Rcpp::_, 0);
        ac += hs1 (Rcpp::_, 1);
    }

    Rcpp::NumericMatrix result (size, 2);
    result (Rcpp::_, 0) = z / (double) ntests;
    result (Rcpp::_, 1) = ac / (double) ntests;
    Rcpp::colnames (result) = Rcpp::CharacterVector::create ("y", "ac");

    return result;
}
