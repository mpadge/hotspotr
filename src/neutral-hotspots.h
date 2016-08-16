#include <Rcpp.h>

#ifndef _RcppNEUTRAL_HOTSPOTS_H
#define _RcppNEUTRAL_HOTSPOTS_H

Rcpp::NumericVector rcpp_trunc_ndist (int len, double sd0);
Rcpp::NumericMatrix rcpp_neutral_hotspots (Rcpp::List nbs, Rcpp::List wts,
        Rcpp::List nbsi, double alpha, double sd0, bool log_scale, int niters,
        std::string ac_type);
Rcpp::NumericMatrix rcpp_neutral_hotspots_ntests (Rcpp::List nbs, 
        Rcpp::List wts, Rcpp::List nbsi, double alpha, double sd0, int niters,
        std::string ac_type, bool log_scale, int ntests);

#endif
