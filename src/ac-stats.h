#include <Rcpp.h>

#ifndef _RcppAC_STATS_H
#define _RcppAC_STATS_H

Rcpp::NumericVector rcpp_ac_stats (Rcpp::NumericVector x, Rcpp::List nbs, 
        Rcpp::List wts, std::string ac_type);

#endif
