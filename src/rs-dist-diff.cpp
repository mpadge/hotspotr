#include <Rcpp.h>
#include "ac-stats.h"
#include "neutral-hotspots.h"

//' rcpp_rs_dist_diff
//'
//' Calculates expected values for deviation of observed rank--scale
//' distributions.
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
//' @param z_mn mean rank-scale distribution of z-variable
//' @param ac_mean mean rank-scale distribution of autocorrelation statistic of
//' z
//'
//' @return A matrix of dimension (size, 2), with first column containing
//' sorted and re-scaled hotspot values, and second column containing sorted and
//' re-scaled spatial autocorrelation statistics.
//'
//' @note \code{rcpp_neutral_hotspots_ntests} returns two vectors containing
//' mean values of raw and autocorrelation statistics from a series of neutral 
//' simulations. This function calculates the distribution of mean squared 
//' differences between these mean profiles and an additional series of simulated
//' instances.
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_rs_dist_diff (Rcpp::List nbs, 
        Rcpp::List wts, double alpha_t, double alpha_s, double sd0, int nt, int
        ntests, std::string ac_type, Rcpp::NumericVector z_mn,
        Rcpp::NumericVector ac_mn)
{
    const int size = nbs.size ();

    double xmin, xmax, tempd;

    Rcpp::NumericVector z (size), ac (size);
    Rcpp::NumericMatrix result (ntests, 2);
    std::fill (result.begin (), result.end (), 0.0);

    for (int n=0; n<ntests; n++)
    {
        z = rcpp_neutral_hotspots (nbs, wts, alpha_t, alpha_s, sd0, nt);
        ac = rcpp_ac_stats (z, nbs, wts, ac_type); // sorted and normalised
        std::sort (z.begin (), z.end (), std::greater<double> ());
        z = (z - (double) Rcpp::min (z)) /
            ((double) Rcpp::max (z) - (double) Rcpp::min (z));

        result (n, 0) = Rcpp::sum ((z - z_mn) * (z - z_mn));
        result (n, 1) = Rcpp::sum ((ac - ac_mn) * (ac - ac_mn));
    }
    Rcpp::colnames (result) = Rcpp::CharacterVector::create ("z", "ac");

    return result;
}
