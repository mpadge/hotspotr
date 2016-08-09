#include <Rcpp.h>
#include <fstream>
#include <algorithm> // std::max

//' rcpp_ives
//'
//' Implements neutral model of Ives & Klopfer (Ecology 1997).
//'
//' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
//' point. 
//' @param alpha_t Strength of temporal autocorrelation
//' @param alpha_s Strength of spatial autocorrelation
//' @param nt Number of successive layers of temporal and spatial autocorrelation
//' used to generate final modelled values
//' @param svec Vector of random numbers for values of \code{s} drawn from a
//' truncated normal ' distribution.
//' @param rvec Vector of random numbers for values of \code{r} drawn from a
//' truncated normal ' distribution.
//'
//' @return A vector of simulated values of same size as \code{nbs}.
//'
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_ives (Rcpp::List nbs, int nt, 
        double alpha_t, double alpha_s,
        Rcpp::NumericVector svec, Rcpp::NumericVector rvec)
{
    const double a0 = 0.05;
    const int size = nbs.size ();

    double tempd;

    if (svec.size () != (size * nt))
        Rcpp::stop ("svec must have length (size * size * nt)");
    if (rvec.size () != (size * nt))
        Rcpp::stop ("rvec must have length (size * size * nt)");

    Rcpp::NumericVector x (size), x2 (size);
    std::fill (x.begin (), x.end (), 1.0);

    Rcpp::NumericVector one_list;
    for (int t=0; t<nt; t++)
    {
        // time step first
        for (int i=0; i<size; i++)
            x (i) = svec [t * size + i] * x (i) * 
                (1.0 + rvec [t * size + i] / (1.0 + a0 * x (i)));
        // spatial autocorrelation
        x2 = Rcpp::clone (x);
        for (int i=0; i<size; i++)
        {
            tempd = 0.0;
            one_list = Rcpp::as <Rcpp::NumericVector> (nbs (i));
            for (int j=0; j<one_list.size (); j++)
                tempd += x (one_list (j) - 1); // nbs are R 1-indexed
            x2 (i) = (1.0 - (double) one_list.size () * alpha_s) * x (i) +
                alpha_s * tempd;
        }
        x = Rcpp::clone (x2);
    }

    return x;
}

//' rcpp_ives_spatial
//'
//' Implements neutral model of Ives & Klopfer (Ecology 1997) with additional
//' spatial ' structure, implented here through replacing generic local
//' autocorrelation ' with movement along maximal local gradients. 
//'
//' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
//' point. 
//' @param alpha_t Strength of temporal autocorrelation
//' @param alpha_s Strength of spatial autocorrelation
//' @param nt Number of successive layers of temporal and spatial autocorrelation
//' used to generate final modelled values
//' @param svec Vector of random numbers for values of \code{s} drawn from a
//' truncated normal ' distribution.
//' @param rvec Vector of random numbers for values of \code{r} drawn from a
//' truncated normal ' distribution.
//'
//' @return A vector of simulated values of same size as \code{nbs}.
//'
//' @section Note: This is not yet implemented! DO NOT USE!
//'
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_ives_spatial (Rcpp::List nbs, int nt, 
        double alpha_t, double alpha_s,
        Rcpp::NumericVector svec, Rcpp::NumericVector rvec)
{
    const double a0 = 0.05;
    const int size = nbs.size ();

    double tempd;

    if (svec.size () != (size * nt))
        Rcpp::stop ("svec must have length (size * size * nt)");
    if (rvec.size () != (size * nt))
        Rcpp::stop ("rvec must have length (size * size * nt)");

    Rcpp::NumericVector x (size), x2 (size);
    std::fill (x.begin (), x.end (), 1.0);

    Rcpp::NumericVector one_list;
    for (int t=0; t<nt; t++)
    {
        // time step first
        for (int i=0; i<size; i++)
            x (i) = svec [t * size + i] * x (i) * 
                (1.0 + rvec [t * size + i] / (1.0 + a0 * x (i)));
        // spatial autocorrelation
        x2 = Rcpp::clone (x);
        for (int i=0; i<size; i++)
        {
            one_list = Rcpp::as <Rcpp::NumericVector> (nbs (i));
            tempd = Rcpp::max (one_list);
            if (tempd < x (i))
                tempd = Rcpp::min (one_list);
            x2 (i) = (1.0 - alpha_s) * x (i) + alpha_s * tempd;
        }
        x = Rcpp::clone (x2);
    }

    return x;
}
