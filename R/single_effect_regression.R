#' @rdname single_effect_regression
#'
#' @title Fitting Single Effect Regression Models
#'
#' @description Methods for fitting single effect regression (SER)
#'   models from summary data.
#'
#' @param nvar An integer no less than 2 specifying the number of
#'   variables.
#'
#' @param An object, possibly a list, specifying the SER prior. The
#'   exact form of this input will depend on which SER model is
#'   used. See \dQuote{details} for further information.
#'
#' @param An object, possibly a list, specifying settings that
#'   affect the model fitting. See \dQuote{details} for further
#'   information.
#'
#' @return An object of class \dQuote{ser} storing key posterior
#'   statistics for all variables and the state of the SER model
#'   fit. It is a list with the following elements:
#'
#' \item{prior}{TO DO: Describe "prior" output here.}
#'
#' \item{posterior}{Describe "posterior" output here.}
#'
#' \item{lbf}{Describe "lbf" output here.}
#'
#' \item{progress}{Describe "progress" output here.}
#' 
#' \item{local}{Describe "local" output here.}
#'
#' \item{global}{Describe "global" output here.}
#'
#' \item{control}{Describe "control" output here.}
#' 
#' @details TO DO: Define what is meant by an "SER model", define
#'   the specific SER models that are implemented, and explain
#'   what are the parameters specific to each model.
#'
#' @note While this interface is technically user-facing, it is
#'   rather basic interface with very few checks for correctness; it
#'   is mainly used to implement the susie model fitting
#'   methods. Use at your own risk.
#'
#' @seealso \code{\link{susie}}
#' 
#' @importFrom utils modifyList
#'
#' @export
#' 
init_ser_model_normal <- function (nvar, 
                                   prior = list(prob = 1/nvar, 
                                                prior_variance = 1),
                                   control = list()) {

  # Basic checks and basic processing of the inputs.
  if (length(prior$prob) == 1)
    prior$prob <- rep(1,nvar)
  prior$prob <- prior$prob/sum(prior$prob)
  control <- modifyList(ser_normal_control_default(),
                          control,keep.null = TRUE)

  # Initialize the single effect regression model.
  fit <- list(prior     = prior,
              lbf       = as.numeric(NA),
              posterior = as.numeric(NA),
              progress  = as.numeric(NA),
              local     = control$estimate_prior_variance &
                          !control$global_prior_variance,
              global    = control$estimate_prior_variance & 
                          !control$global_prior_variance,
              control   = control)
  class(fit) <- c("ser_normal","ser","list")
  return(fit)
}

#' @rdname single_effect_regression
#' 
#' @export
#' 
ser_normal_control_default <- function()
  list(estimate_prior_variance = TRUE,
       global_prior_variance = FALSE,
       paramtol = 1e-6,
       lbftol = 1e-6,
       maxiter = 10,
       version = "R",
       verbose = "none")

#' @rdname single_effect_regression
#'
#' @param fit An object of class \dQuote{ser} that captures the state
#'   of an SER model fit, typically the output from running an
#'   \dQuote{init_ser_model} function, or an output from an
#'   \code{fit_ser_model} call.
#' 
#' @param dat A matrix with \code{nvar} rows and 2 named columns:
#'   \dQuote{bhat}, the least-squares estimates of the coefficients;
#'   \dQuote{s}, the variances of these estimates. See
#'   \dQuote{details} for how these are precisely defined.
#'
#' @param return_data If \code{return_data = TRUE}, the input data are
#'   also included in the output.
#' 
#' @export
#'
fit_ser_model <- function (fit, dat, return_data = FALSE) {
   UseMethod("fit_ser_model")
}

#' @rdname single_effect_regression
#'
#' @export
#' 
fit_ser_model.ser_normal <- function (fit, dat, return_data = FALSE) {
  # TO DO: 
  # - Write a function to simulate different X matrices---
  #   see https://github.com/stephenslab/dsc-linreg.
  # - Add a "progress" output.
  # - Test the ability to "refit".
  # - Implement an Rcpp version of this.
  # - Add some simple unit tests.
  # - Benchmark the R and C++ implementations.
  # - Fill out the docs, and in particular define the SER model,
  #   and give a small example illustrating its use.
  # - Run R CMD check ...
  if (missing(dat)) {
    if (is.null(fit$dat))
      stop("data not provided",domain = NA)
  }

  # Estimate the prior variance parameter locally.
  if (fit$control$local & fit$control$estimate_prior_variance) {
    for (iter in 1:fit$control$maxiter) {
      prior_variance_current <- fit$prior$prior_variance
      lbf_current <- fit$lbf

      # E step.    
      res <- compute_posterior_ebnm_normal(dat[,"bhat"],dat[,"s"],
                                           fit$prior$prior_variance)
      ans <- normalizelogweights(log(fit$prior$prob) + res$lbf)

      # M-step.
      fit$prior$prior_variance <- sum(ans$probs * (res$b1^2 + res$s1))
      if (fit$control$verbose == "detailed")
      cat("s0 =",format(fit$prior$prior_variance,digits = 12),
          " lbf =",format(ans$weight_sum,digits = 12),"\n")
      if (!is.na(lbf_current) && lbf_current > ans$weight_sum + 1e-8)
        stop("prior variance update decreased the likelihood",domain = NA)
      if (abs(prior_variance_current - 
              fit$prior$prior_variance) < fit$control$paramtol)
        break
    }
  }

  # Computing the posterior statistics one last time and prepare the
  # output.
  res <- compute_posterior_ebnm_normal(dat[,"bhat"],dat[,"s"],
                                       fit$prior$prior_variance)
  ans <- normalizelogweights(log(fit$prior$prob) + res$lbf)
  pip <- ans$probs
  post_mean <- pip * res$b1
  fit$posterior <- cbind(lbf       = res$lbf,
                         pip       = pip,
                         post_mean = post_mean)
  fit$lbf <- ans$weight_sum 
  if (return_data)
    fit$dat <- dat
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
