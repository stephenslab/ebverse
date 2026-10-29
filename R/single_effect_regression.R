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
#' @param dat An object of class \sQuote{ser_data}, typically an
#'     output from \code{prepare_ser_data}.
#'
#' @param prior An object, possibly a list, specifying the SER
#'   prior. The exact form of this input will depend on which SER
#'   model is used. See \sQuote{details} for further information.
#'
#' @param control An object, possibly a list, specifying settings that
#'   affect the model fitting. See \sQuote{details} for further
#'   information.
#'
#' @param return_data If \code{return_data = TRUE}, the input data are
#'   included in the output.
#' 
#' @return \code{init_ser_model} and \code{fit_ser_model} return
#'   objects of class \sQuote{ser} storing key posterior statistics
#'   for all variables and the state of the SER model fit. It is a
#'   list with the following elements:
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
#' \item{progress}{A data frame containing detailed information about
#'   the algorithm's progress.}
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
#' @examples
#' # TO DO: Add an example.
#' 
#' @export
#' 
init_ser_model <- function (dat, 
                            family = "normal",
                            prior = list(prob = 1/nrow(dat), 
                                         prior_variance = 1),
                            control = list(),
                            return_data = FALSE) {
  family <- match.arg(family)
  if (family == "normal")
    fit <- init_ser_model_normal(dat,prior,control,return_data)
  return(fit)
}

#' @rdname single_effect_regression
#' 
#' @param bhat A vector of least-squares estimates of the coefficients.
#'
#' @param shat A vector of variances of the least-squares estimates.
#'
#' @param X A matrix of with rows corresponding to observations
#'   (samples) and columns corresponding to candidate explanatory
#'   variables.
#'
#' @param y A vector of regression outcomes.
#' 
#' @param sigma Residual variances.
#' 
#' @return \code{prepare_ser_data} returns a matrix with two named
#'   columns: \sQuote{bhat}, the least-squares estimates of the
#'   coefficients; \sQuote{shat}, the variances of these
#'   estimates. See \sQuote{details} for precise definitions of
#'   these terms.
#'   
#' @export
#' 
prepare_ser_data <- function (bhat, shat, X, y, sigma) {
  if (!((missing(X) & missing(y) & missing(sigma) & 
         !missing(bhat) & !missing(shat)) |
        (!missing(X) & !missing(y) & !missing(sigma) & 
         missing(bhat) & missing(shat))))
    stop("Either {bhat, shat} or {X, y, sigma} should be provided ",
         "(and not both)")
  if (missing(bhat) & missing(shat)) {
    xx   <- colSums(X^2)
    bhat <- drop(y %*% X)/xx
    shat <- sigma/xx
  }
  if (length(bhat) != length(shat))
    stop("bhat and s should be the same length")
  if (length(bhat) < 2 | length(shat) < 2)
    stop("bhat and s should have at least two entries")
  dat <- cbind(bhat = bhat,shat = shat)
  class(dat) <- c("ser_data","matrix","array")
  return(dat)
}

#' @rdname single_effect_regression
#'
#' @importFrom utils modifyList
#'
#' @export
#' 
init_ser_model_normal <- function (dat, 
                                   prior = list(prob = 1/nrow(dat), 
                                                prior_variance = 1),
                                   control = list(),
                                   return_data = FALSE) {

  # Basic checks and processing of the inputs.
  if (!inherits(dat,"ser_data"))
    stop("dat is not a valid SER data set")
  nvar <- nrow(dat)
  if (length(prior$prob) != 1 & length(prior$prob) != nvar)
    stop("prior$prob should be a vector of length 1 or nvar",domain = NA)
  if (length(prior$prob) == 1)
    prior$prob <- rep(1,nvar)
  if (prior$prior_variance <= 0)
    stop("prior$prior_variance should be a positive number",domain = NA)
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
  if (return_data)
    fit$dat <- dat
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
       lbftol = as.numeric(NA),
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
#' @export
#'
fit_ser_model <- function (fit, dat, return_data = missing(dat)) {
   UseMethod("fit_ser_model")
}

#' @rdname single_effect_regression
#'
#' @export
#' 
fit_ser_model.ser_normal <- function (fit, dat, return_data = missing(dat)) {

  # TO DO: 
  #
  # - Write a function to simulate different X matrices---
  #   see https://github.com/stephenslab/dsc-linreg.
  # - Implement an Rcpp version of this.
  # - Add some simple unit tests.
  # - Benchmark the R and C++ implementations.
  # - Fill out the docs, and in particular define the SER model,
  #   and give a small example illustrating its use.
  # - Run R CMD check...
  #
  force(return_data)
  if (missing(dat)) {
    if (is.null(fit$dat))
      stop("the data were not provided",domain = NA)
    dat <- fit$dat
  }  

  # Estimate the prior variance parameter locally.
  res <- fit_ser_model_normal_core(dat,fit$prior,fit$lbf,fit$control)
  fit$prior     <- res$prior
  fit$lbf       <- res$lbf
  fit$posterior <- res$posterior
  progress      <- res$progress

  # Prepare the output.
  if (is.data.frame(fit$progress)) {
    progress$iter <- progress$iter + max(fit$progress$iter)
    fit$progress <- rbind(fit$progress,progress)
  } else
    fit$progress <- res$progress
  if (return_data)
    fit$dat <- dat
  else
    fit["dat"] <- NULL
  return(fit)
}

# This implements the core model fitting algorithm for the single
# effect regression with normal priors.
fit_ser_model_normal_core <- function (dat, prior, lbf, control) {
  progress <- as.numeric(NA)
  if (control$local & control$estimate_prior_variance) {
    progress <- as.data.frame(matrix(0,control$maxiter,5))
    names(progress) <- c("iter","lbf","prior.var","delta.lbf",
                         "delta.prior.var")
    progress[,"iter"] <- 1:control$maxiter
    if (control$verbose == "detailed")
      cat("iter log-Bayes factor prior var d(logBF) d(p.var)\n")
    for (iter in 1:control$maxiter) {
      prior_variance_current <- prior$prior_variance
      lbf_current <- lbf

      # E step.    
      res <- compute_posterior_ebnm_normal(dat[,"bhat"],dat[,"shat"],
                                           prior$prior_variance)
      ans <- normalizelogweights(log(prior$prob) + res$lbf)
      lbf <- ans$weight_sum

      # M-step.
      prior$prior_variance <- sum(ans$probs * (res$b1^2 + res$s1))

      # Update progress data frame and check stopping criteria.
      progress[iter,"lbf"] <- lbf
      progress[iter,"prior.var"] <- prior$prior_variance
      progress[iter,"delta.lbf"] <- lbf - lbf_current
      progress[iter,"delta.prior.var"] <- 
        prior$prior_variance - prior_variance_current
      if (control$verbose == "detailed")
        cat(sprintf("%4d %+16.9e %8.3e %+8.1e %+8.1e\n",
                    progress[iter,"iter"],progress[iter,"lbf"],
                    progress[iter,"prior.var"],progress[iter,"delta.lbf"],
                    progress[iter,"delta.prior.var"]))
      if (prior$prior_variance <= 0)
        stop("prior variance is zero or negative",domain = NA)
      if (!is.na(lbf_current) && 
          !is.na(control$lbftol) &&
          lbf_current > lbf + control$lbftol)
        stop("prior variance update decreased the likelihood",domain = NA)
      if (!is.na(control$paramtol) &&
          abs(prior_variance_current - 
              prior$prior_variance) < control$paramtol)
        break
      if (!is.na(lbf_current) && 
          !is.na(control$lbftol) &&
          lbf - lbf_current < control$lbftol)
        break
    }
  }

  # Compute the posterior statistics one last time and prepare the
  # output.
  res <- compute_posterior_ebnm_normal(dat[,"bhat"],dat[,"shat"],
                                       prior$prior_variance)
  ans <- normalizelogweights(log(prior$prob) + res$lbf)
  pip <- ans$probs
  b2  <- res$s1 + res$b1^2
  post_mean <- pip*res$b1
  post_var  <- pip*b2 - post_mean^2
  posterior <- cbind(lbf  = res$lbf,
                     pip  = pip,
                     mean = post_mean,
                     var  = post_var)
  lbf       <- ans$weight_sum 
  progress  <- progress[1:iter,]
  rownames(progress) <- NULL
  return(list(prior = prior,lbf = lbf,posterior = posterior,
              progress = progress))
}

# Compute the Bayes factor and posterior statistics for an EBNM model
# with a zero-centered normal prior.
#
#' @importFrom stats dnorm
compute_posterior_ebnm_normal <- function (bhat, s, s0) {
  lbf <- dnorm(bhat,0,sqrt(s0 + s),log = TRUE) - 
         dnorm(bhat,0,sqrt(s),log = TRUE)
  s1 <- 1/(1/s0 + 1/s)
  b1 <- bhat*s1/s
  return(list(lbf = lbf,b1 = b1,s1 = s1))
}
