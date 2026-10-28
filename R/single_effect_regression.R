#' @rdname single_effect_regression
#'
#' @title Fitting Single Effect Regression Models
#'
#' @description Methods for fitting single effect regression (SER)
#'   models and computing posterior statistics from summary
#'   data. Typically the model fitting involves estimating
#'   parameters in the prior (\dQuote{hyperparameters}) by
#'   maximizing the likelihood, but potentially other (non-prior)
#'   parameters could be estimated as well. The \sQuote{fit_ser_model}
#'   method is the main method for fitting SER models, and for
#'   computing posterior statistics.
#'
#' @param nvar An integer no less than 2 specifying the number of
#'   variables.
#'
#' @param prior An object, possibly a list, specifying the SER
#'   prior. The exact form of this input will depend on which SER
#'   model is used. See \sQuote{details} for further information.
#'
#' @param control An object, possibly a list, specifying settings that
#'   affect the model fitting. See \sQuote{details} for further
#'   information.
#'
#' @return An object of class \sQuote{ser} storing key posterior
#'   statistics for all variables and the state of the SER model
#'   fit. It is a list with the following elements:
#'
#' \item{prior}{Data structure specifying the prior. This may
#'   include hyperparameters that change during model fitting.}
#'
#' \item{posterior}{A matrix with one row per variable and the
#'   following columns: \sQuote{lbf}, the log-Bayes factor;
#'   \sQuote{PIP}, the posterior inclusion probability; \sQuote{mean},
#'   the posterior mean of the regression coefficient; \sQuote{var},
#'   the posterior variance of the regression coefficient.}
#'
#' \item{lbf}{Log-Bayes factor comparing the SER model to the
#'   \dQuote{null} model (b = 0).}
#'
#' \item{progress}{TO DO: Describe "progress" output here.}
#' 
#' \item{local}{Logical specifying whether one or more parameters may
#'   be estimated \dQuote{locally} (that is, they are estimated
#'   separately in this SER).}
#'
#' \item{global}{Logical specifying whether one or more parameters may
#'   be estimated \dQuote{globally} (that is, they are estimated
#'   simultaneously in multiple SERs).}
#'
#' \item{control}{A data structure, likely a list, containing the
#'   settings that affect how the model fitting is conducted.}
#' 
#' @details TO DO: Define what is meant by an "SER model", define
#'   the specific SER models that are implemented, and explain
#'   what are the parameters specific to each model.
#'
#' @note This is a basic interface with very few checks; it is mainly
#'   used to implement the susie model fitting methods. Use at your
#'   own risk.
#'
#' @seealso \code{\link{fit_susie_model}}
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
  if (length(prior$prob) != 1 & length(prior$prob) != nvar)
    stop("prior$prob should be a vector of length 1 or nvar")
  if (length(prior$prob) == 1)
    prior$prob <- rep(1,nvar)
  if (prior$prior_variance <= 0)
    stop("prior$prior_variance should be a positive number")
  prior$prob <- prior$prob/sum(prior$prob)
  control <- modifyList(ser_normal_control_default(),
                          control,keep.null = TRUE)

  # Initialize the single effect regression model.
  fit <- list(prior     = prior,
              lbf       = as.numeric(NA),
              posterior = as.numeric(NA),
              progress  = as.numeric(NA),
              control   = control)
  fit$control$local <- control$estimate_prior_variance &
                       !control$global_prior_variance
  fit$control$global <- control$estimate_prior_variance & 
                        !control$global_prior_variance
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
#' @param fit An object of class \sQuote{ser} that captures the state
#'   of an SER model fit, typically the output from running an
#'   \sQuote{init_ser_model} function, or an output from a
#'   \code{fit_ser_model} call.
#' 
#' @param dat A matrix with \code{nvar} rows and 2 named columns:
#'   \sQuote{bhat}, the least-squares estimates of the coefficients;
#'   \sQuote{s}, the variances of these estimates. See
#'   \sQuote{details} for how these are precisely defined.
#'
#' @param return_data If \code{return_data = TRUE}, the input data are
#'   included in the output.
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
  # - Implement tolerance checking.
  # - Test the ability to "refit".
  # - Implement an Rcpp version of this.
  # - Add some simple unit tests.
  # - Benchmark the R and C++ implementations.
  # - Fill out the docs, and in particular define the SER model,
  #   and give a small example illustrating its use.
  if (missing(dat)) {
    if (is.null(fit$dat))
      stop("the data were not provided",domain = NA)
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
  b2  <- res$s1 + res$b1^2
  post_mean <- pip*res$b1
  post_var <-  pip*b2 - post_mean^2
  fit$posterior <- cbind(lbf  = res$lbf,
                         pip  = pip,
                         mean = post_mean,
                         var  = post_var)
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
