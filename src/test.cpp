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
  Rcpp::IntegerVector extract_vector(int i);
  int get_genotype(int i, int j);
private:
  pcaMatrix(const pcaMatrix&);
  pcaMatrix& operator=(const pcaMatrix&);
  boost::interprocess::file_mapping file;
  boost::interprocess::mapped_region file_region;
  const char* file_data;
  const int nrow;
  const int ncol;
};

pcaMatrix::pcaMatrix(std::string path, const int n, const int p) : nrow(n), ncol(p) {
  try {
    this->file = file_mapping(path.c_str(), read_only);
  } catch(const interprocess_exception& e) {
    throw std::runtime_error("File not found.");
  }
  this->file_region = mapped_region(this->file, read_only);
  this->file_data = static_cast<const char*>(this->file_region.get_address());
}

int pcaMatrix::get_genotype(int i, int j) {
  // Reduce two-dimensional index to one-dimensional index with the mode
  int which_pos = (j * this->nrow) + i;
  return which_pos;
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
RcppExport int pcaMatrix__get_genotype(SEXP xp_, SEXP i_, SEXP j_) {
  // Convert inputs to appropriate C++ types
  Rcpp::XPtr<pcaMatrix> ptr(xp_);
  int i = Rcpp::as<int>(i_);
  int j = Rcpp::as<int>(j_);
  try {
    // Invoke the extract_vector function
    int res = ptr->get_genotype(i, j);
    return res;
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
    return 0;
  }
};