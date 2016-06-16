#include <Rcpp.h>
#include <fstream>

// [[Rcpp::export]]
Rcpp::NumericMatrix rcpp_ives2D (int size, int nt, 
        double alpha_t, double alpha_s,
        Rcpp::NumericVector svec, Rcpp::NumericVector rvec)
{
    int indx;
    const double s0 = 0.5, r0 = 1.05, a0 = 0.05;
    double s, r;

    Rcpp::NumericMatrix x (size, size), x2 (size, size);
    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++)
            x (i, j) = 1.0;

    //std::ofstream out_file;
    //out_file.open ("000.txt", std::ofstream::out);
    //out_file.close ();

    for (int t=0; t<nt; t++)
    {
        // time step first
        for (int i=0; i<size; i++)
            for (int j=0; j<size; j++)
            {
                indx = t * size * size + i * size + j;
                x (i, j) = svec [indx] * x (i, j) * 
                    (1.0 + rvec [indx] / (1.0 + a0 * x(i, j)));
            }
        x2 = Rcpp::clone (x);
        // then spatial autocorrelation
        for (int i=1; i<(size - 1); i++)
            for (int j=1; j<(size - 1); j++)
                x2 (i, j) = (1.0 - 4 * alpha_t) * x (i ,j) + alpha_t * 
                    (x(i-1, j) + x(i+1, j) + x(i, j-1) + x(i, j+1));
        // wrap edges
        for (int i=1; i<(size - 1); i++)
        {
            x2 (i, 0) = (1.0 - 4 * alpha_s) * x (i ,0) + alpha_s * 
                (x(i-1, 0) + x(i+1, 0) + x(i, size-1) + x(i, 1));
            x2 (i, size-1) = (1.0 - 4 * alpha_s) * x (i ,size-1) + alpha_s * 
                (x(i-1, size-1) + x(i+1, size-1) + x(i, size-2) + x(i, 0));
            x2 (0, i) = (1.0 - 4 * alpha_s) * x (0 ,i) + alpha_s * 
                (x(0, i-1) + x(0, i+1) + x(size-1, i) + x(1, i));
            x2 (size-1, i) = (1.0 - 4 * alpha_s) * x (size-1, i) + alpha_s * 
                (x(size-1, i-1) + x(size-1, i+1) + x(size-2, i) + x(0, i));
        }
        // and corners
        x2 (0, 0) = (1.0 - 4 * alpha_s) * x (0, 0) + alpha_s *
            (x (0, 1) + x (0, size-1) + x (size-1, 0) + x (1, 0));
        x2 (0, size-1) = (1.0 - 4 * alpha_s) * x (0, size-1) + alpha_s *
            (x (0, size-2) + x (0, 0) + x (1, size-1) + x (size-1, size-1));
        x2 (size-1, 0) = (1.0 - 4 * alpha_s) * x (size-1, 0) + alpha_s *
            (x (size-1, 1) + x (size-1, size-1) + x (0, 0) + x (size-2, 0));
        x2 (size-1, size-1) = (1.0 - 4 * alpha_s) * x (size-1, size-1) + alpha_s *
            (x (size-2, size-1) + x (0, size-1) + x (size-1, 0) + x (size-1, size-2));

        x = Rcpp::clone (x2);
    }

    return x;
}
