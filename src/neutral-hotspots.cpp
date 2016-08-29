#include <Rcpp.h>
#include "ac-stats.h"
#include <fstream>

//' rcpp_trunc_ndist
//'
//' Truncated normal distribution (mean 1, respective upper and lower limits of
//' 0 and 2). Code copied directly from `github.com/mpadge/tnorm`, with the
//' readme of that repo demonstrating the speed advantages of using this rather
//' than pre-existing approaches (the R package `truncnorm`).
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

    std::vector <double> z;
    z.resize (0);
    while (z.size () < len)
    {
        eps = Rcpp::rnorm (len, 1.0, sd);
        for (Rcpp::NumericVector::iterator it = eps.begin ();
                it != eps.end (); ++it)
            if (*it >= 0.0 && *it <= 2.0)
                z.push_back (*it);
    }
    z.resize (len);

    std::copy (z.begin (), z.end (), eps.begin ());
    z.resize (0);

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
    Rcpp::NumericVector z1 = rcpp_trunc_ndist (size, sd0); 
    Rcpp::NumericVector z2 (size), nbs_to, nbs_from, nbs_n, ac;

    // Spatial autocorrelation
    for (int i=0; i<niters; i++)
    {
        std::fill (z2.begin (), z2.end (), 0.0);
        for (int j=0; j<nbsi.size (); j++)
        {
            tempmat = Rcpp::as <Rcpp::NumericMatrix> (nbsi (j));
            nbs_to = tempmat (Rcpp::_, 0);
            nbs_from = tempmat (Rcpp::_, 1);
            nbs_n = tempmat (Rcpp::_, 2);

            // Note that nbs are 1-indexed
            for (int k=0; k<nbs_to.size (); k++)
            {
                z2 (nbs_to (k) - 1) += ((1.0 - alpha) * z1 (nbs_to (k) - 1) +
                    alpha * z1 (nbs_from (k) - 1)) / (double) nbs_n (k);
            }
        } // end for j over nbsi
        std::copy (z2.begin (), z2.end (), z1.begin ());
    } // end for i over niters

    if (log_scale)
    {
        // negative values sometimes arise for alpha < 0, but scale is
        // arbitrary, so re-scaled here to ensure all values are > 0
        z1 = 1.0 + z1 - (double) Rcpp::min (z1);
        z1 = log10 (z1);
    }

    ac = rcpp_ac_stats (z1, nbs, wts, ac_type);
    std::sort (z1.begin (), z1.end (), std::greater<double> ());
    
    z1 = (z1 - (double) Rcpp::min (z1)) /
        ((double) Rcpp::max (z1) - (double) Rcpp::min (z1));

    Rcpp::NumericMatrix result (size, 2);
    result (Rcpp::_, 0) = z1;
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
    Rcpp::colnames (result) = Rcpp::CharacterVector::create ("z", "ac");

    return result;
}
