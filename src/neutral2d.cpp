#include <Rcpp.h>

const double DOUBLE_MAX = std::numeric_limits<double>::max (),
    DOUBLE_MIN = -DOUBLE_MAX;

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_neutral2d_1test (int size, 
        double alpha_t, double alpha_s, double sd0, int nt)
{
    int indx;
    double tempd;

    // Set up truncated normal distribution
    Rcpp::NumericVector eps (nt * size * size);
    eps = Rcpp::rnorm (nt * size * size, 1.0, sd0);
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
        

    Rcpp::NumericMatrix x (size, size), x2 (size, size);
    std::fill (x.begin (), x.end (), 1.0);

    for (int t=0; t<nt; t++)
    {
        // temporal autocorrelation
        for (int i=1; i<size; i++)
            for (int j=0; j<size; j++)
            {
                indx = t * size * size + i * size + j;
                x (i, j) = (1.0 - alpha_t) * x (i, j) + alpha_t * eps (indx);
            }
        x2 = Rcpp::clone (x);
        // spatial autocorrelation
        for (int i=1; i<(size - 1); i++)
            for (int j=1; j<(size - 1); j++)
                x2 (i, j) = (1.0 - 4.0 * alpha_s) * x (i, j) + alpha_s *
                    (x (i-1, j) + x (i+1, j) + x (i, j-1) + x (i, j+1));
        // edges
        for (int i=1; i<(size - 1); i++)
        {
            x2 (i, 0) = (1.0 - 4.0 * alpha_s) * x (i, 0) + alpha_s *
                (x (i-1, 0) + x (i+1, 0) + x (i, size-1) + x (i, 1));
            x2 (i, size-1) = (1.0 - 4.0 * alpha_s) * x (i, size-1) + alpha_s *
                (x (i-1, size-1) + x (i+1, size-1) + x (i, size-2) + x (i, 0));
            x2 (0, i) = (1.0 - 4.0 * alpha_s) * x (0, i) + alpha_s *
                (x (size-1, i) + x (1, i) + x (0, i-1) + x (0, i+1));
            x2 (size-1, i) = (1.0 - 4.0 * alpha_s) * x (size-1, i) + alpha_s *
                (x (size-1, i-1) + x (size-1, i+1) + x (size-2, i) + x (0, i));
        }
        // corners
        x2 (0, 0) = (1.0 - 4.0 * alpha_s) * x (0, 0) + alpha_s *
            (x (0, 1) + x (0, size-1) + x (size-1, 0) + x (1, 0));
        x2 (0, size-1) = (1.0 - 4.0 * alpha_s) * x (0, size-1) + alpha_s *
            (x (0, size-2) + x (0, 0) + x (1, size-1) + x (size-1, size-1));
        x2 (size-1, 0) = (1.0 - 4.0 * alpha_s) * x (size-1, 0) + alpha_s *
            (x (size-1, 1) + x (size-1, size-1) + x (0, 0) + x (size-2, 0));
        x2 (size-1, size-1) = (1.0 - 4.0 * alpha_s) * x (size-1, size-1) + 
                alpha_s * (x (size-2, size-1) + x (0, size-1) + 
                x (size-1, size-2) + x (size-1, 0));

        x = Rcpp::clone (x2);
    } // end for t

    return x;
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_morani (Rcpp::NumericMatrix x)
{
    int size = x.nrow ();
    Rcpp::NumericMatrix xexp (size + 2, size + 2);
    // Can't figure out how to use Rcpp::Range to set, rather than extract, a
    // sub-matrix, so nmat is filled explicitly here
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
            tempd = xexp (i-1, j) + xexp (i+1, j) +
                xexp (i, j-1) + xexp (i, j+1);
            mmat (i-1, j-1) = (xexp (i, j) - xmn) * (tempd - 4.0 * xmn) / 5.0;
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

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_neutral2d_ntests (int size, 
        double alpha_t, double alpha_s, double sd0, int nt, int ntests)
{
    int len = size * size;
    double xmin, xmax;

    Rcpp::NumericMatrix x (size, size);

    Rcpp::NumericVector xv (len), xv_tot (len), ac_vec (len), ac_vec1 (len);
    std::fill (ac_vec.begin (), ac_vec.end (), 0.0);
    std::fill (xv_tot.begin (), xv_tot.end (), 0.0);
    for (int n=0; n<ntests; n++)
    {
        x = rcpp_neutral2d_1test (size, alpha_t, alpha_s, sd0, nt);
        ac_vec1 = rcpp_morani (x); // already sorted and normalised
        for (int i=0; i<size; i++)
            for (int j=0; j<size; j++)
                xv (i * size + j) = x (i, j);
        for (int i=0; i<len; i++)
            ac_vec (i) = ac_vec (i) + ac_vec1 (i);
        std::sort (xv.begin (), xv.end (), std::greater<double> ());
        xmin = Rcpp::min (xv);
        xmax = Rcpp::max (xv);
        for (int i=0; i<len; i++)
            xv_tot (i) = xv_tot (i) + (xv (i) - xmin) / (xmax - xmin);
    }

    Rcpp::NumericMatrix result (len, 2);
    for (int i=0; i<len; i++)
    {
        result (i, 0) = xv_tot (i) / (double) ntests;
        result (i, 1) = ac_vec (i) / (double) ntests;
    }

    return result;
}
