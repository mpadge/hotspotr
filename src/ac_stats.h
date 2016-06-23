#include <Rcpp.h>

#ifndef _RcppAC_STATS_H
#define _RcppAC_STATS_H

Rcpp::NumericVector rcpp_ac_stats (Rcpp::List nbs, Rcpp::NumericVector x, 
        std::string ac_type);

#endif
