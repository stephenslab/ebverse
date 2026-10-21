#include "RcppArmadillo.h"

// funcPtr will take a const reference to a vec and return a vector.
typedef arma::vec (*funcPtr)(const arma::vec& x);

