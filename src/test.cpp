// [[Rcpp::depends(BH)]]

#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/exceptions.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <Rcpp.h>
#include <string>

using namespace boost::interprocess;

class pcaMatrix {
public:
  pcaMatrix(std::string path, int n, int p);
  int get_genotype(std::size_t i, std::size_t j);
  Rcpp::IntegerVector extract_vector(Rcpp::IntegerVector i);
  Rcpp::IntegerMatrix extract_matrix(Rcpp::IntegerVector i, Rcpp::IntegerVector j);
  Rcpp::NumericVector prodvect(Rcpp::NumericVector x, int nIND, int nSNP);
  Rcpp::NumericVector crossprodvect(Rcpp::NumericVector x, int nIND, int nSNP);
private:
  pcaMatrix(const pcaMatrix&);
  pcaMatrix& operator=(const pcaMatrix&);
  boost::interprocess::file_mapping file;
  boost::interprocess::mapped_region file_region;
  const char* file_data;
  int nrow;
  int ncol;
};

pcaMatrix::pcaMatrix(std::string path, int n, int p) : nrow(n), ncol(p) {
  try {
    this->file = file_mapping(path.c_str(), read_only);
  } catch(const interprocess_exception& e) {
    throw std::runtime_error("File not found.");
  }
  this->file_region = mapped_region(this->file, read_only);
  this->file_data = static_cast<const char*>(this->file_region.get_address());
}

int pcaMatrix::get_genotype(std::size_t i, std::size_t j) {
  std::size_t which_pos = (j * this->nrow) + i;
  char genotype = this->file_data[2 * which_pos];
  int mapping = NA_INTEGER; // missing
  if (genotype == 48) {
    mapping = 0; // homozygous AA
  } else if (genotype == 49) {
    mapping = 1; // homozygous BB
  } else if (genotype == 50) {
    mapping = 2; // heterozygous AB
  }
  return mapping;
}

Rcpp::IntegerVector pcaMatrix::extract_vector(Rcpp::IntegerVector i) {
  // Convert from 1-index to 0-index
  Rcpp::IntegerVector i0(i - 1);
  // Keep size of i
  std::size_t size_i = i.size();
  // Reserve output vector
  Rcpp::IntegerVector out(size_i);
  // Get bounds
  std::size_t bounds = this->nrow * this->ncol;
  // Iterate over indexes
  for (std::size_t idx_i = 0; idx_i < size_i; idx_i++) {
    if (Rcpp::IntegerVector::is_na(i0[idx_i]) || static_cast<std::size_t>(i0[idx_i]) >= bounds) {
      out(idx_i) = NA_INTEGER;
    } else {
      out(idx_i) = this->get_genotype(i0[idx_i] % this->nrow, i0[idx_i] / this->nrow);
    }
  }
  return out;
}

Rcpp::IntegerMatrix pcaMatrix::extract_matrix(Rcpp::IntegerVector i, Rcpp::IntegerVector j) {
  // Check if indexes are out of bounds
  if (Rcpp::is_true(Rcpp::any(i > this->nrow)) || Rcpp::is_true(Rcpp::any(j > this->ncol))) {
    throw std::runtime_error("subscript out of bounds");
  }
  // Convert from 1-index to 0-index
  Rcpp::IntegerVector i0(i - 1);
  Rcpp::IntegerVector j0(j - 1);
  // Keep sizes of i and j
  std::size_t size_i = i.size();
  std::size_t size_j = j.size();
  // Reserve output matrix
  Rcpp::IntegerMatrix out(size_i, size_j);
  // Iterate over column indexes
  for (std::size_t idx_j = 0; idx_j < size_j; idx_j++) {
    for (std::size_t idx_i = 0; idx_i < size_i; idx_i++) {
      out(idx_i, idx_j) = this->get_genotype(i0[idx_i], j0[idx_j]);
    }
  }
  return out;
}

Rcpp::NumericVector pcaMatrix::prodvect(Rcpp::NumericVector x, int nIND, int nSNP){
  Rcpp::NumericVector res(nIND);
  // Iterate over column indexes
  for (int i = 0; i < nIND; i++) {
    for (int j = 0; j < nSNP; j++) {
      res[i] += (double) (x[j] * this-> get_genotype(i, j));
    }
  }
  return(res);
} 

Rcpp::NumericVector pcaMatrix::crossprodvect(Rcpp::NumericVector x, int nIND, int nSNP){
  Rcpp::NumericVector res(nSNP);
  // Iterate over column indexes
  for (int j = 0; j < nSNP; j++) {
    for (int i = 0; i < nIND; i++) {
      res[j] += (double) (x[i] * this-> get_genotype(i, j));
    }
  }
  return(res);
}

// Export BEDMatrix::BEDMatrix
RcppExport SEXP pcaMatrix__new(SEXP path_, SEXP n_, SEXP p_) {
  // Convert inputs to appropriate C++ types
  std::string path = Rcpp::as<std::string>(path_);
  int n = Rcpp::as<int>(n_);
  int p = Rcpp::as<int>(p_);
  try {
    // Create a pointer to a BEDMatrix object and wrap it as an external
    // pointer
    Rcpp::XPtr<pcaMatrix> ptr(new pcaMatrix(path, n, p), true);
    // Return the external pointer to the R side
    return ptr;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
    return 0;
  }
};

// Export BEDMatrix::extract_vector
RcppExport SEXP pcaMatrix__extract_vector(SEXP xp_, SEXP i_) {
  // Convert inputs to appropriate C++ types
  Rcpp::XPtr<pcaMatrix> ptr(xp_);
  Rcpp::IntegerVector i = Rcpp::as<Rcpp::IntegerVector>(i_);
  try {
    // Invoke the extract_vector function
    Rcpp::IntegerVector res = ptr->extract_vector(i);
    return res;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
    return 0;
  }
};

// Export BEDMatrix::extract_matrix
RcppExport SEXP pcaMatrix__extract_matrix(SEXP xp_, SEXP i_, SEXP j_) {
  // Convert inputs to appropriate C++ types
  Rcpp::XPtr<pcaMatrix> ptr(xp_);
  Rcpp::IntegerVector i = Rcpp::as<Rcpp::IntegerVector>(i_);
  Rcpp::IntegerVector j = Rcpp::as<Rcpp::IntegerVector>(j_);
  try {
    // Invoke the extract_matrix function
    Rcpp::IntegerMatrix res = ptr->extract_matrix(i, j);
    return res;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
    return 0;
  }
};

// Export BEDMatrix::extract_matrix
RcppExport SEXP pcaMatrix__prodvect(SEXP xp_, SEXP x_, SEXP nIND_, SEXP nSNP_) {
  // Convert inputs to appropriate C++ types
  Rcpp::XPtr<pcaMatrix> ptr(xp_);
  Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(x_);
  int nIND = Rcpp::as<int>(nIND_);
  int nSNP = Rcpp::as<int>(nSNP_);
  try {
    // Invoke the extract_matrix function
    Rcpp::NumericVector res = ptr->prodvect(x, nIND, nSNP);
    return res;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
    return 0;
  }
};

// Export BEDMatrix::extract_matrix
RcppExport SEXP pcaMatrix__crossprodvect(SEXP xp_, SEXP x_, SEXP nIND_, SEXP nSNP_) {
  // Convert inputs to appropriate C++ types
  Rcpp::XPtr<pcaMatrix> ptr(xp_);
  Rcpp::NumericVector x = Rcpp::as<Rcpp::NumericVector>(x_);
  int nIND = Rcpp::as<int>(nIND_);
  int nSNP = Rcpp::as<int>(nSNP_);
  try {
    // Invoke the extract_matrix function
    Rcpp::NumericVector res = ptr->crossprodvect(x, nIND, nSNP);
    return res;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
    return 0;
  }
};




