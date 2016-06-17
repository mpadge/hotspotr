#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_neutral1d (int size, double alpha_t, double alpha_s,
        int nt, Rcpp::NumericVector eps)
{
    int indx;

    size = size * size;
    Rcpp::NumericVector y (size);
    for (int i=0; i<size; i++)
        y (i) = 1.0;

    for (int t=0; t<nt; t++)
    {
        // temporal autocorrelation
        for (int i=0; i<size; i++)
        {
            indx = t * size + i;
            y (i) = (1.0 - alpha_t) * y (i) + alpha_t * eps (indx);
        }
        // spatial autocorrelation
        for (int i=1; i<(size - 1); i++)
            y (i) = (1.0 - 2.0 * alpha_s) * y (i) +
                    alpha_s * (y (i-1) + y (i+1));
    }

    return y;
}
