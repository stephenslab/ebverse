#include "ebverse_types.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List rcpp_hello_world() {
  CharacterVector x = CharacterVector::create("foo","bar");
  NumericVector   y = NumericVector::create(0.0,1.0);
  List            z = List::create(x,y);
  return z;
}

// Source: gallery.rcpp.org/articles/passing-cpp-function-pointers/
arma::vec fun1_cpp(const arma::vec& x) {
  arma::vec y = x + x;
  return y;
}

arma::vec fun2_cpp(const arma::vec& x) {
  arma::vec y = 10*x;
  return y;
}

//
// The following function takes a string argument, picks a function
// and returns it wrapped as an external pointer SEXP. We could return
// this to R as well.
//
// [[Rcpp::export]]
XPtr<funcPtr> putFunPtrInXPtr(std::string fstr) {
  if (fstr == "fun1")
    return(XPtr<funcPtr>(new funcPtr(&fun1_cpp)));
  else if (fstr == "fun2")
    return(XPtr<funcPtr>(new funcPtr(&fun2_cpp)));
  else
    return XPtr<funcPtr>(R_NilValue);
}

// [[Rcpp::export]]
arma::vec callViaXPtr(const arma::vec x, SEXP xpsexp) {
  XPtr<funcPtr> xpfun(xpsexp);
  funcPtr fun = *xpfun;
  arma::vec y = fun(x);
  return (y);
}
