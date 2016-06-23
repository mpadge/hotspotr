#include <Rcpp.h>
#include "ac_stats.h"

//' rcpp_ac_stats
//'
//' Computes spatial autocorrelation statistics for a given input matrix
//'
//' @param nbs An \code{spdep} \code{nb} object listing all neighbours of each
//' point
//' @param x Corresponding vector of values
//' @param ac_type Character string specifying type of aucorrelation
//' (\code{moran}, \code{geary}, or code{getis-ord}).
//'
//' @return A vector of sorted spatial autocorrelation statistics scaled between
//' zero and one.
//'
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_ac_stats (Rcpp::List nbs, Rcpp::NumericVector x, 
        std::string ac_type)
{
    if (x.size () != nbs.size ())
        Rcpp::stop ("nbs must be same size as x");

    const int size = x.size ();
    const double xmn = Rcpp::mean (x), sd=Rcpp::sd (x);

    double tempd [2], sizei;

    Rcpp::NumericVector ac (size), one_list;

    for (int i=0; i<size; i++)
    {
        one_list = Rcpp::as <Rcpp::NumericVector> (nbs (i));
        // nbs are 1-indexed from R, so converted here to 0-indexed
        for (auto j = one_list.begin (); j != one_list.end (); ++j)
            (*j) = (*j) - 1;

        tempd [0] = tempd [1] = 0.0;
        if (ac_type == "moran")
        {
            for (auto j = one_list.begin (); j != one_list.end (); ++j)
                tempd [0] += (x (*j) - xmn) * (x (*j) - xmn);
            for (auto j = one_list.begin (); j != one_list.end (); ++j)
                for (auto k = j; ++k != one_list.end (); /**/)
                    tempd [1] += (x (*j) - xmn) * (x (*k) - xmn);
            ac (i) = tempd [1] / tempd [0];
        } else if (ac_type == "geary")
        {
            for (auto j = one_list.begin (); j != one_list.end (); ++j)
                tempd [0] += (x (*j) - xmn) * (x (*j) - xmn);
            for (auto j = one_list.begin (); j != one_list.end (); ++j)
                for (auto k = j; ++k != one_list.end (); /**/)
                    tempd [1] += (x (*j) - x (*k)) * (x (*j) - x (*k));
            ac (i) = tempd [1] / tempd [0];
        } else // getis-ord
        {
            for (auto j = one_list.begin (); j != one_list.end (); ++j)
                tempd [0] += x (*j);
            sizei = (double) one_list.size ();
            tempd [0] = tempd [0] - sizei * xmn;
            tempd [1] = sizei * ((double) size - sizei);
            tempd [1] = sqrt (tempd [1] / ((double) size - 1.0));
            //ac (i) = tempd [0] / (sd * tempd [1]);
            ac (i) = tempd [0] / tempd [1];
        }
    }

    std::sort (ac.begin (), ac.end (), std::greater<double> ());
    double amax = Rcpp::max (ac), amin = Rcpp::min (ac);
    for (Rcpp::NumericVector::iterator it = ac.begin (); it != ac.end (); ++it)
        *it = (*it - amin) / (amax - amin);

    return ac;
}
