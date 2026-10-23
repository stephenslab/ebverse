#' @rdname single_effect_regression
#'
#' @title Fitting Single Effect Regression Models
#'
#' @description Add description here.
#'
#' @param nvar Describe the "nvar" input here.
#'
#' @param Describe the "prior_family" input here.
#' 
#' @param Describe the "prior" input here.
#'
#' @param Describe the "control" input here.
#'
#' @return An object storing: key posterior statistics for all
#'     variables; the state of the model fit. It is a list
#'     with the following elements:
#'
#' \item{prior}{Describe "prior" output here.}
#'
#' \item{posterior}{Describe "posterior" output here.}
#'
#' \item{control}{Describe "control" output here.}
#'
#' @details Explain that this is a "minimal" user interface with very
#'     little input checking; use directly at your own risk. It is
#'     mainly used to implement the susie methods.
#' 
#' @export
#' 
init_ser_model <- function (nvar, prior_family, prior, control) {
  # TO DO.
}

#' @importFrom utils modifyList
#'
#' @export
#' 
init_ser_model_normal <- function (nvar, 
                                   prior = list(prob = 1/nvar, 
                                                prior_variance = 1),
                                   control = list()) {
    control <- modifyList(ser_normal_control_default(),
                        control,keep.null = TRUE)
  if (length(prior$prob) == 1)
    prior$prob <- rep(1,nvar)
  prior$prob <- prior$prob/sum(prior$prob)
  fit <- list(prior = prior,
              posterior = as.numeric(NA),
              progress = as.numeric(NA),
              control = control)
  class(fit) <- c("ser_normal","ser","list")
  return(fit)
}

#' @rdname single_effect_regression
#' 
#' @export
#' 
ser_normal_control_default <- function()
  list(estimate_prior_variance = TRUE,
       tol = 1e-6,
       maxiter = 10,
       version = "R",
       verbose = "none")

#' @rdname single_effect_regression
#'
#' @export
#'
fit_ser_model <- function (fit, dat, return_data = FALSE) {
   UseMethod("fit_ser_model")
}

#' @rdname single_effect_regression
#'
#' @param dat Describe the "dat" input here.
#' 
#' @export
#' 
fit_ser_model.ser_normal <- function (fit, dat, return_data = FALSE) {
  # TO DO: 
  # - Implement return_data option.
  # - Implement estimation of the prior variance.
  # - Create a "simulate" S3 method.
  # - Write a function to simulate different X matrices---
  #   see https://github.com/stephenslab/dsc-linreg.
  # - Add a "progress" output.
  # - Implement "verbose" option.
  # - Test the ability to "refit".
  # - Implement an Rcpp version of this.
  # - Add some simple unit tests.
  # - Benchmark the R and C++ implementations.
  # - Implement Laplace prior ("ser_laplace").
  # - Fill out the docs, and in particular define the SER model,
  #   and give a small example illustrating its use.
  w <- fit$prior$prob
  res <- compute_posterior_ebnm_normal(dat[,"bhat"],dat[,"s"],
                                       fit$prior$prior_variance)
  ans <- normalizelogweights(log(w) + res$lbf)
  pip <- ans$probs
  posterior_mean <- pip * res$b1
  fit$posterior <- cbind(lbf = res$lbf,
                         pip = pip,
                         posterior_mean = posterior_mean)
  fit$lbf <- ans$weight_sum
  return(fit)
}

# TO DO: Explain here what this function is for, and what are the
# inputs and outputs.
#
#' @importFrom stats dnorm
compute_posterior_ebnm_normal <- function (bhat, s, s0) {
  lbf <- dnorm(bhat,0,sqrt(s0 + s),log = TRUE) - 
         dnorm(bhat,0,sqrt(s),log = TRUE)
  s1 <- 1/(1/s0 + 1/s)
  b1 <- bhat*s1/s
  return(list(lbf = lbf,b1 = b1,s1 = s1))
}
