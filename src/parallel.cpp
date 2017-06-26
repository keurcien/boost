#include <Rcpp.h>
#include <omp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector sum_vec_1(NumericVector x, NumericVector y) {
  NumericVector z(x.size());
  #pragma omp parallel for
  for (int i = 0; i < x.size(); i++){
    z[i] = x[i] + y[i];
  }
  return z;
}


// [[Rcpp::export]]
NumericVector sum_vec_3(NumericVector x, NumericVector y) {
  NumericVector z(x.size());
  for (int i = 0; i < x.size(); i++){
    z[i] = x[i] + y[i];
  }
  return z;
}


