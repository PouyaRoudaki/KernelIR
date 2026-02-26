#ifndef HANDLE_DUPLICATES_AND_TIES_H
#define HANDLE_DUPLICATES_AND_TIES_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Function declaration (prototype)
Rcpp::IntegerMatrix handle_duplicates_and_ties(
    Rcpp::NumericMatrix X,
    Rcpp::IntegerMatrix nn_idx,
    Rcpp::NumericMatrix nn_dist,
    int Knn
);

#endif // HANDLE_DUPLICATES_AND_TIES
