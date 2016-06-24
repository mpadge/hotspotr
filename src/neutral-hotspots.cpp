#include <Rcpp.h>
#include "ac-stats.h"

//' rcpp_neutral_hotspots
//'
//' Implements neutral model in two dimensions
//'
//' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
//' point. 
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
Rcpp::NumericVector rcpp_neutral_hotspots (Rcpp::List nbs, 
        double alpha_t, double alpha_s, double sd0, int nt)
{
    const int size = nbs.size ();

    int indx;
    double tempd;

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
        

    Rcpp::NumericVector x (size), x2 (size);
    std::fill (x.begin (), x.end (), 1.0);

    Rcpp::NumericVector one_list;
    for (int t=0; t<nt; t++)
    {
        // time step first
        for (int i=0; i<size; i++)
            x (i) = (1.0 - alpha_t) * x (i) + alpha_t * eps (t * size + i);
        // spatial autocorrelation
        x2 = Rcpp::clone (x);
        for (int i=0; i<size; i++)
        {
            tempd = 0.0;
            one_list = Rcpp::as <Rcpp::NumericVector> (nbs (i));
            for (int j=0; j<one_list.size (); j++)
                tempd += x (j);
            x2 (i) = (1.0 - (double) one_list.size () * alpha_s) * x (i) +
                alpha_s * tempd;
        }
        x = Rcpp::clone (x2);
    }

    return x;
}


//' rcpp_neutral_hotspots_ntests
//'
//' Performs repeated neutral tests to yield average distributions of both
//' hotspot values and spatial autocorrelation statistics.
//'
//' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
//' point. 
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
        double alpha_t, double alpha_s, double sd0, int nt, int ntests,
        std::string ac_type)
{
    const int size = nbs.size ();

    Rcpp::NumericVector x (size), xtot (size), ac (size), ac1 (size);
    std::fill (ac.begin (), ac.end (), 0.0);
    std::fill (xtot.begin (), xtot.end (), 0.0);

    for (int n=0; n<ntests; n++)
    {
        x = rcpp_neutral_hotspots (nbs, alpha_t, alpha_s, sd0, nt);
        ac1 = rcpp_ac_stats (nbs, x, ac_type); // sorted and normalised
        for (int i=0; i<size; i++)
            ac (i) += ac1 (i);
        std::sort (x.begin (), x.end (), std::greater<double> ());
        for (int i=0; i<size; i++)
            xtot (i) += (x (i) - (double) Rcpp::min (x)) /
                ((double) Rcpp::max (x) - (double) Rcpp::min (x));
    }

    Rcpp::NumericMatrix result (size, 2);
    for (int i=0; i<size; i++)
    {
        result (i, 0) = xtot (i) / (double) ntests;
        result (i, 1) = ac (i) / (double) ntests;
    }
    Rcpp::colnames (result) = Rcpp::CharacterVector::create ("y", "ac");

    return result;
}
