#include <Rcpp.h>

#ifndef _RcppNEUTRAL_HOTSPOTS_H
#define _RcppNEUTRAL_HOTSPOTS_H

Rcpp::NumericVector rcpp_neutral_hotspots (Rcpp::List nbs, Rcpp::List wts,
        double alpha_t, double alpha_s, double sd0, int nt);

#endif
