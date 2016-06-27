#include <Rcpp.h>

#ifndef _RcppNEUTRAL_HOTSPOTS_H
#define _RcppNEUTRAL_HOTSPOTS_H

Rcpp::NumericVector rcpp_neutral_hotspots (Rcpp::List nbs, Rcpp::List wts,
        double alpha_t, double alpha_s, double sd0, int nt);
Rcpp::NumericMatrix rcpp_neutral_hotspots_ntests (Rcpp::List nbs, 
        Rcpp::List wts, double alpha_t, double alpha_s, double sd0, int nt, int
        ntests, std::string ac_type);

#endif
