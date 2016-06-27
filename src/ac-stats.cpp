#include <Rcpp.h>
#include "ac-stats.h"

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
Rcpp::NumericVector rcpp_ac_stats (Rcpp::NumericVector z, Rcpp::List nbs, 
        Rcpp::List wts, std::string ac_type)
{
    if (z.size () != nbs.size ())
        Rcpp::stop ("nbs must be same size as z");

    const int size = z.size ();
    const double zmn = Rcpp::mean (z), sd=Rcpp::sd (z);

    double tempd [4], sizei;

    Rcpp::NumericVector ac (size), nbs1, wts1;

    for (int i=0; i<size; i++)
    {
        nbs1 = Rcpp::as <Rcpp::NumericVector> (nbs (i));
        wts1 = Rcpp::as <Rcpp::NumericVector> (wts (i));
        // nbs are 1-indexed from R, so converted here to 0-indexed
        for (auto j = nbs1.begin (); j != nbs1.end (); ++j)
            (*j) = (*j) - 1;

        tempd [0] = tempd [1] = tempd [2] = tempd [3] = 0.0;
        if (ac_type == "moran")
        {
            for (int j=0; j<nbs1.size (); j++)
            {
                tempd [0] += (z (nbs1 (j)) - zmn) * (z (nbs1 (j)) - zmn) *
                    wts1 (j);
                tempd [2] += wts1 (j);
            }
            for (int j=0; j<(nbs1.size () - 1); j++)
                for (int k=(j+1); k<nbs1.size (); k++)
                {
                    tempd [1] += (z (nbs1 (j)) - zmn) * (z (nbs1 (k)) - zmn) *
                        wts1 (j) * wts1 (k);
                    tempd [3] += wts1 (j) * wts1 (k);
                }
            ac (i) = tempd [1] * tempd [2] / (tempd [0] * tempd [3]);
        } else if (ac_type == "geary")
        {
            for (int j=0; j<nbs1.size (); j++)
            {
                tempd [0] += (z (nbs1 (j)) - zmn) * (z (nbs1 (j)) - zmn) *
                    wts1 (j);
                tempd [2] += wts1 (j);
            }
            for (int j=0; j<(nbs1.size () - 1); j++)
                for (int k=(j+1); k<nbs1.size (); k++)
                {
                    tempd [1] += (z (nbs1 (j)) - nbs1 (k)) * 
                        (z (nbs1 (j)) - nbs1 (k)) * wts1 (j) * wts1 (k);
                    tempd [3] += wts1 (j) * wts1 (k);
                }
            ac (i) = tempd [1] * tempd [2] / (tempd [0] * tempd [3]);
        } else // getis-ord
        {
            for (int j=0; j<nbs1.size (); j++)
            {
                tempd [0] += z (nbs1 (j)) * wts1 (j);
                tempd [2] += wts1 (j);
            }
            sizei = (double) nbs1.size ();
            tempd [0] = tempd [0] / tempd [2] - sizei * zmn;
            tempd [1] = sizei * ((double) size - sizei);
            tempd [1] = sqrt (tempd [1] / ((double) size - 1.0));
            //ac (i) = tempd [0] / (sd * tempd [1]);
            ac (i) = tempd [0] / tempd [1];
        }
    }

    std::sort (ac.begin (), ac.end (), std::greater<double> ());
    ac = (ac - Rcpp::min (ac)) / (Rcpp::max (ac) - Rcpp::min (ac));

    return ac;
}
