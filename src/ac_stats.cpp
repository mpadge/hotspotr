#include <Rcpp.h>
#include "ac_stats.h"

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_ac_stats (Rcpp::NumericMatrix x, std::string ac_type)
{
    int size = x.nrow ();
    Rcpp::NumericMatrix xexp (size + 2, size + 2);
    // Can't figure out how to use Rcpp::Range to set, rather than extract, a
    // sub-matrix, so nmat is filled explicitly here
    // 
    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
            xexp (i + 1, j + 1) = x (i, j);
    xexp (Rcpp::_, 0) = xexp (Rcpp::_, size);
    xexp (0, Rcpp::_) = xexp (size, Rcpp::_);
    xexp (Rcpp::_, size+1) = xexp (Rcpp::_, 1);
    xexp (size+1, Rcpp::_) = xexp (1, Rcpp::_);

    Rcpp::NumericMatrix mmat (size, size);
    double tempd, xmn = mean (x);
    for (int i=1; i<(size+1); i++)
        for (int j=1; j<(size+1); j++)
        {
            if (ac_type == "moran")
            {
                tempd = xexp (i-1, j) + xexp (i+1, j) +
                    xexp (i, j-1) + xexp (i, j+1);
                mmat (i-1, j-1) = (xexp (i, j) - xmn) * (tempd - 4.0 * xmn) / 5.0;
            } else if (ac_type == "geary")
            {
                mmat (i-1, j-1) = 
                    (xexp (i, j) - xexp (i-1, j)) * (xexp (i, j) - xexp (i-1, j)) +
                    (xexp (i, j) - xexp (i+1, j)) * (xexp (i, j) - xexp (i+1, j)) +
                    (xexp (i, j) - xexp (i, j-1)) * (xexp (i, j) - xexp (i, j-1)) +
                    (xexp (i, j) - xexp (i, j+1)) * (xexp (i, j) - xexp (i, j+1));
            } else
            {
                mmat (i-1, j-1) = xexp (i-1, j) + xexp (i+1, j) + xexp (i, j-1) + 
                                    xexp (i, j+1) - 4.0 * xmn;
            }
        }
    Rcpp::NumericVector mvec (size * size);
    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
            mvec (i * size + j) = mmat (i, j);
    std::sort (mvec.begin (), mvec.end (), std::greater<double> ());
    double mmin = Rcpp::min (mvec);
    double mmax = Rcpp::max (mvec);
    for (Rcpp::NumericVector::iterator it = mvec.begin (); it != mvec.end (); ++it)
        *it = (*it - mmin) / (mmax - mmin);

    return mvec;
}
