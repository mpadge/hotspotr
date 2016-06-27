#include <Rcpp.h>
#include "ac-stats.h"
#include "neutral-hotspots.h"

//' rcpp_p_values
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
//'
//' @return A matrix of dimension (size, 2), with first column containing
//' sorted and re-scaled hotspot values, and second column containing sorted and
//' re-scaled spatial autocorrelation statistics.
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_p_values (Rcpp::List nbs, 
        Rcpp::List wts, double alpha_t, double alpha_s, double sd0, int nt, int
        ntests, std::string ac_type)
{
    const int size = nbs.size ();

    double xmin, xmax, tempd;

    // Generate initial mean values from 100 tests
    Rcpp::NumericMatrix test_means = rcpp_neutral_hotspots_ntests (nbs, wts,
            alpha_t, alpha_s, sd0, nt, 100, ac_type);

    Rcpp::NumericVector x (size), ac (size);
    Rcpp::NumericMatrix result (ntests, 4);
    std::fill (result.begin (), result.end (), 0.0);

    for (int n=0; n<ntests; n++)
    {
        x = rcpp_neutral_hotspots (nbs, wts, alpha_t, alpha_s, sd0, nt);
        ac = rcpp_ac_stats (nbs, wts, x, ac_type); // sorted and normalised
        std::sort (x.begin (), x.end (), std::greater<double> ());
        xmin = (double) Rcpp::min (x);
        xmax = (double) Rcpp::max (x);
        for (int i=0; i<size; i++)
            x (i) = (x (i) - xmin) / (xmax - xmin);

        for (int i=0; i<size; i++)
        {
            tempd = x (i) - test_means (i, 0);
            result (n, 0) += fabs (tempd);
            result (n, 1) += tempd * tempd;
            tempd = ac (i) - test_means (i, 1);
            result (n, 2) += fabs (tempd);
            result (n, 3) += tempd * tempd;
        }
        result (n, 0) = result (n, 0) / (double) size;
        result (n, 1) = result (n, 1) / (double) size - 
            result (n, 0) * result (n, 0);
        result (n, 1) = result (n, 1) * (double) size / ((double) size - 1.0);
        result (n, 2) = result (n, 2) / (double) size;
        result (n, 3) = result (n, 3) / (double) size - 
            result (n, 2) * result (n, 2);
        result (n, 3) = result (n, 3) * (double) size / ((double) size - 1.0);
    }
    Rcpp::colnames (result) = Rcpp::CharacterVector::create ("z", "z2",
            "ac", "ac2");

    return result;
}
