#' @name single_effect_regression
#' 
#' @title Fitting Single Effect Regression Models
#'
#' @description Methods for fitting single effect regression (SER)
#'   models and computing posterior statistics from summary
#'   data. Typically the model fitting involves estimating
#'   parameters in the prior (\dQuote{hyperparameters}) by
#'   maximizing the likelihood, but potentially other parameters
#'   could be estimated as well. The \sQuote{fit_ser_model} method
#'   is the main method for fitting SER models and returning key
#'   posterior statistics. SER models are implemented as a modular
#'   (\dQuote{object-oriented}) interface, allowing for different
#'   types of SER families (see \sQuote{Details} for specification
#'   of SER classes).
#'
#' @param variable.names The names of the variables in the
#'   SER model.
#' 
#' @param family The SER family. See \dQuote{Details} for more
#'   information.
#' 
#' @param prior An object, possibly a list, specifying the SER
#'   prior. The exact form of this input will depend on which SER
#'   model is used. See \sQuote{Details} for further information.
#'
#' @param init_method Posterior statistics need to be in SER model
#'   fits at all times, even if they are just initial guesses. This
#'   argument specifies how the initial guess is obtained.
#' 
#' @param control An object, possibly a list, specifying settings that
#'   affect the model fitting. See \sQuote{Details} for further
#'   information.
#'
#' @return \code{init_ser_model} and \code{fit_ser_model} return
#'   objects of class \sQuote{ser} storing key posterior statistics
#'   for all variables and the state of the SER model fit. It is a
#'   list with the following elements:
#'
#'   \item{prior}{Data structure specifying the prior. This may
#'     include hyperparameters that change during model fitting.}
#'
#'   \item{posterior}{A data frame with one row per variable and the
#'     following columns: \sQuote{lbf}, the log-Bayes factor;
#'     \sQuote{PIP}, the posterior inclusion probability; \sQuote{mean},
#'     the posterior mean of the regression coefficient; \sQuote{var},
#'     the posterior variance of the regression coefficient.}
#'
#'   \item{lbf}{Log-Bayes factor comparing the SER model to the
#'     \dQuote{null} model (b = 0).}
#'
#'   \item{KL}{Kullback-Leibler divergence between the SER posterior
#'     distribution and the SER prior.}
#'
#'   \item{progress}{A data frame containing detailed information about
#'     the progress on the model fitting.}
#' 
#'   \item{local}{Logical specifying whether one or more parameters may
#'     be estimated \dQuote{locally} (that is, they are estimated
#'     separately in this SER).}
#'
#'   \item{global}{Logical specifying whether one or more parameters
#'     may be estimated \dQuote{globally} (that is, they are estimated
#'     simultaneously in multiple SERs). This is used when the SER
#'     model is used as a building block for another model such as a
#'     SuSiE model.}
#'
#'   \item{control}{A data structure, likely a list, containing the
#'     settings that affect how the model fitting is conducted.}
#'
#'   \item{dat}{When \code{return_data = TRUE}, this contains the 
#'     data that were used to fit the model.}
#' 
#' @details A single effect regression (SER) model is a Bayesian
#'   multiple regression model in which only a single candidate
#'   variable is active (\dQuote{included}). Therefore, the
#'   posterior distribution becomes a \dQuote{mixture} of simpler
#'   single-variable regression models.
#'
#'   The first element of an SER model is a multiple regression
#'   model of an outcome y,
#'   \deqn{y = Xb + e,}
#'   where \eqn{e} is a vector of independent, normally
#'   distributed error terms, \eqn{e \sim N(0, \sigma^2 I)}.
#'   The residual variance \eqn{\sigma^2} is assumed to be known.
#'
#'   The second element of a an SER model is an \emph{SER prior} on
#'   the coefficients vector \eqn{b}. An SER prior has the
#'   property that any \eqn{b} drawn from the prior has exactly
#'   one nonzero coefficient. For example, the SER model initialized
#'   by \code{init_ser_model_normal()} defines a prior in which
#'   coefficient \eqn{j} is nonzero with prior probability
#'   \eqn{\pi_j}, and the nonzero coefficient is assigned a normal
#'   prior with zero mean and a variance \eqn{\sigma_0^2} which may be
#'   estimated. Because the SER posterior reduces to a mixture of
#'   single-variable regression models, the sufficient statistics
#'   needed to compute the posterior are the same as the sufficient
#'   statistics needed to compute the posterior under each
#'   single-variable model, such as the least-squares estimate of each
#'   coefficient and its variance.
#' 
#'   The key steps of an analysis with a single effect regression are:
#'   (i) prepare the data using the \code{prepare_ser_data} function;
#'   (ii) initialize the SER model using the \code{init_ser_data}
#'   function; and (iii) fit the model and compute posterior
#'   statistics from the fitted model using the \code{fit_ser_data}
#'   function. 
#'
#'   \bold{Design of single effect regression classes.} A
#'   \dQuote{single effect regression} (SER) class has at a minimum
#'   the following three elements:
#'   \enumerate{
#' 
#'   \item A function that initializes the SER model (e.g.,
#'     \code{init_ser_model_normal}). It returns an \sQuote{ser}
#'     model.
#' 
#'   \item A \sQuote{fit_ser_model} method that takes (at least) two
#'     inputs, an \sQuote{ser} model and an \sQuote{ser_data} data
#'     set. It returns an updated SER model (object of class
#'     \sQuote{ser}).
#'
#'   \item A \sQuote{variable.names} method that returns the
#'     variable names. 
#'   }
#' 
#' @note This is a lightweight interface with relatively few features
#'   and few checks; it is mainly used to implement the susie
#'   methods.
#'
#' @seealso \link{susie}
#'
#' @examples
#' # TO DO: Add one or several examples.
#' 
#' @export
#' 
init_ser_model <- function (variable.names, 
                            family = "normal",
                            prior = list(prob = 1/length(variable.names),
                                         prior_variance = 1),
                            init_method = "null",
                            control = list()) {
  # TO DO: Implement "init_method" option for posterior means to be
  # specified manually.
  family <- match.arg(family)
  if (family == "normal")
    fit <- init_ser_model_normal(variable.names,prior,init_method,control)
  return(fit)
}

#' @rdname single_effect_regression
#'
#' @importFrom utils modifyList
#' @importFrom stats rnorm
#' 
#' @export
#' 
init_ser_model_normal <- 
  function (variable.names, 
            prior = list(prob = 1/length(variable.names), 
                         prior_variance = 1),
            init_method = c("null","random"),
            control = list()) {

  # Basic checks and processing of the inputs.
  nvar <- length(variable.names)
  if (length(prior$prob) != 1 & length(prior$prob) != nvar)
    stop("prior$prob should be a vector of length 1 or nvar",domain = NA)
  if (length(prior$prob) == 1)
    prior$prob <- rep(1,nvar)
  if (prior$prior_variance <= 0)
    stop("prior$prior_variance should be a positive number",domain = NA)
  prior$prob <- prior$prob/sum(prior$prob)
  names(prior$prob) <- variable.names
  init_method <- match.arg(init_method)
  control <- modifyList(ser_normal_control_default(),control,
                        keep.null = TRUE)

  # Make an initial guess for the posterior statistics.
  posterior <- data.frame(lbf  = 0,
                          pip  = prior$prob,
                          mean = 0,
                          var  = prior$prior_variance)
  rownames(posterior) <- variable.names
  if (init_method == "random")
    posterior$mean <- rnorm(nvar,0,sqrt(prior$prior_variance))

  # Initialize the SER model.
  fit <- list(prior     = prior,
              lbf       = as.numeric(NA),
              KL        = as.numeric(NA),
              posterior = posterior,
              progress  = as.numeric(NA),
              control   = control)
  fit$control$local  <- control$estimate_prior_variance &
                        !control$global_prior_variance
  fit$control$global <- control$estimate_prior_variance & 
                        control$global_prior_variance
  class(fit) <- c("ser_normal","ser","list")
  return(fit)
}

#' @rdname single_effect_regression
#'
#' @param fit An object of class \sQuote{ser} that captures the state
#'   of an SER model fit, typically the output from running an
#'   \sQuote{init_ser_model} function, or an output from a previous
#'   \code{fit_ser_model} call.
#' 
#' @param return_data If \code{return_data = TRUE}, the input data are
#'   included in the output.
#'
#' @param \dots Additional arguments controlling the model fitting.
#' 
#' @export
#'
fit_ser_model <- function (fit, dat, return_data = missing(dat), ...) {
   UseMethod("fit_ser_model")
}

#' @rdname single_effect_regression
#'
#' @param dat An object of class \sQuote{ser_data}, such as the
#'   output from a call to \code{prepare_ser_data}.
#'
#' @export
#' 
fit_ser_model.ser_normal <- function (fit, dat, return_data = missing(dat), 
                                      ...) {

  # TO DO: 
  # - Implement function to conveniently simulate susie data sets
  #   of different sizes.
  # - Run benchmarking on different sizes of J and compare with 
  #   susie(L = 1) and susie_rss(L = 1).
  # - Implement an Rcpp version of this and see if it speeds up 
  #   fit_susie_model() for the small n, large J case.
  # - Add some simple unit tests.
  # - Benchmark the R and C++ implementations.
  # - Fill out the docs, and in particular define the SER model,
  #   and give a small example illustrating its use.
  # - Run R CMD check...
  force(return_data)
  if (missing(dat)) {
    if (is.null(fit$dat))
      stop("the data were not provided",domain = NA)
    dat <- fit$dat
  }  

  # Fit the SER model and compute the posterior statistics.
  res <- fit_ser_model_normal_core(dat,fit$prior,fit$lbf,fit$control)
  fit$prior     <- res$prior
  fit$lbf       <- res$lbf
  fit$KL        <- res$KL
  fit$posterior <- res$posterior
  progress      <- res$progress

  # Prepare the output.
  rownames(fit$posterior) <- rownames(dat)
  if (fit$control$trace) {
    if (is.data.frame(fit$progress)) {
      progress$iter <- progress$iter + max(fit$progress$iter)
      fit$progress <- rbind(fit$progress,progress)
    } else
      fit$progress <- progress
  } else
    fit$progress <- as.numeric(NA)
  if (return_data)
    fit$dat <- dat
  else
    fit["dat"] <- NULL
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
       lbftol   = as.numeric(NA),
       maxiter  = 10,
       version  = "R",
       trace    = TRUE,
       verbose  = "none")

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
#' @param residual_variance The residual variance parameter which is
#'   used to compute the the outputted statistics. This is only used
#'   when X, y are provided instead of bhat, shat.
#'
#' @note In \code{prepare_set_data}, only provide (bhat, shat) or
#'   (X, y); an error will occur if the wrong combination of data
#'   are provided.
#' 
#' @return \code{prepare_ser_data} returns a data frame with two
#'   columns: \sQuote{bhat}, the least-squares estimates of the
#'   coefficients; \sQuote{shat}, the variances of these
#'   estimates. See \sQuote{Details} for precise definitions of
#'   these terms.
#'   
#' @export
#' 
prepare_ser_data <- function (bhat, shat, X, y, residual_variance = 1) {
  if (!((missing(X) & missing(y) & !missing(bhat) & !missing(shat)) |
        (!missing(X) & !missing(y) & missing(bhat) & missing(shat))))
    stop("Either {bhat, shat} or {X, y} should be provided ",
         "(and not both)",domain = NA)
  if (missing(bhat) & missing(shat)) {
    J <- ncol(X)
    if (length(residual_variance) != 1)
      stop("residual_variance should be a single number",domain = NA)
    res  <- compute_lse(X,y,residual_variance)
    bhat <- res$bhat
    shat <- res$shat
  } else
    residual_variance <- as.numeric(NA)
  if (length(bhat) != length(shat))
    stop("bhat and s should be the same length",domain = NA)
  if (length(bhat) < 2 | length(shat) < 2)
    stop("bhat and s should have at least two entries",domain = NA)
  if (is.null(names(bhat)) | is.null(names(shat)))
    stop("variables must be named by providing names for bhat and shat, or ",
         "column names for X",domain = NA)
  dat <- data.frame(bhat = bhat,shat = shat)
  rownames(dat) <- names(bhat)
  class(dat) <- c("ser_data","data.frame")
  attr(dat,"residual_variance") <- residual_variance
  return(dat)
}

# This implements the core model fitting algorithm for the single
# effect regression with normal priors.
fit_ser_model_normal_core <- function (dat, prior, lbf, control) {
  progress <- as.numeric(NA)

  # Estimate the prior variance parameter "locally".
  if (control$estimate_prior_variance & control$local) {
    if (control$trace | control$verbose == "detailed") {
      progress <- as.data.frame(matrix(0,control$maxiter + 1,5))
      names(progress) <- c("iter","lbf","prior.var","d.lbf","d.prior.var")
      progress[,"iter"] <- seq(1,control$maxiter + 1)
      if (control$verbose == "detailed") {
        cat("iter log-Bayes factor   sigma0 d(logBF) d(sigma0)\n")
        cat("---- ---------------- -------- -------- ---------\n")
      }
    }
    for (iter in 1:control$maxiter) {

      # Save some quantities that are used to monitor progress.
      current_state <- list(prior_variance = prior$prior_variance,lbf = lbf)

      # E step.
      post <- compute_posterior_ebnm_normal(dat$bhat,dat$shat,
                                            prior$prior_variance)
      ans <- normalizelogweights(log(prior$prob) + post$lbf)
      lbf <- ans$weight_sum
      pip <- ans$probs
      b1  <- post$b1
      s1  <- post$s1

      # M-step.
      prior$prior_variance <- sum(pip * (b1^2 + s1))
      
      # Update the progress data frame.
      if (is.data.frame(progress)) {
        progress[iter,"lbf"]         <- lbf
        progress[iter,"prior.var"]   <- prior$prior_variance
        progress[iter,"d.lbf"]       <- lbf - current_state$lbf
        progress[iter,"d.prior.var"] <- prior$prior_variance - 
                                        current_state$prior_variance
      }

      # Update the user about the algorithm's process.
      if (control$verbose == "detailed")
        cat(sprintf("%4d %+16.9e %08.1e %+8.1e %+8.2e\n",
                    progress[iter,"iter"],
                    progress[iter,"lbf"],
                    progress[iter,"prior.var"],
                    progress[iter,"d.lbf"],
                    progress[iter,"d.prior.var"]))

      # Check for errors as well as the stopping criteria. Some checks
      # are only performed after the first iteration.
      if (prior$prior_variance <= 0)
        stop("prior variance is zero or negative",domain = NA)
      if (iter > 1) {
        if (!is.na(current_state$lbf) && 
            !is.na(control$lbftol) &&
            lbf < current_state$lbf)
          warning("prior variance update decreased the likelihood",
                  domain = NA)
        if (!is.na(control$paramtol) &&
            abs(current_state$prior_variance - 
                prior$prior_variance) < control$paramtol) {
          if (control$verbose != "none")
            message("parameter tolerance is met; model fitting is stopped",
                    domain = NA)
          break
        }
        if (!is.na(current_state$lbf) && 
            !is.na(control$lbftol) &&
            lbf - current_state$lbf < control$lbftol) {
          if (control$verbose != "none")  
            message("Bayes factor tolerance is met; model fitting is ",
                    "stopped",domain = NA)
          break
        }
      }
    }
  }

  # Compute the posterior statistics one last time.
  current_state <- list(lbf = lbf)
  post <- compute_posterior_ebnm_normal(dat$bhat,dat$shat,
                                        prior$prior_variance)
  ans <- normalizelogweights(log(prior$prob) + post$lbf)
  lbf <- ans$weight_sum 
  pip <- ans$probs
  b1  <- post$b1
  s1  <- post$s1
  b2  <- s1 + b1^2
  post_mean <- pip*b1
  post_var  <- pip*b2 - post_mean^2

  # Compute the K-L divergence between the posterior and prior.
  KL <- compute_kl_ser(dat$bhat,dat$shat,post_mean,post_var,lbf)

  # Update the progress data frame to reflect the new posterior
  # calculations.
  if (control$trace) {
    iter <- iter + 1
    progress[iter,"lbf"]         <- lbf
    progress[iter,"prior.var"]   <- prior$prior_variance
    progress[iter,"d.lbf"]       <- lbf - current_state$lbf
    progress[iter,"d.prior.var"] <- 0
    progress <- progress[1:iter,]
  } else
    progress <- as.numeric(NA)

  # Prepare the output.
  posterior <- data.frame(lbf  = post$lbf,
                          pip  = pip,
                          mean = post_mean,
                          var  = post_var)
  return(list(prior     = prior,
              lbf       = lbf,
              KL        = KL,
              posterior = posterior,
              progress  = progress))
}

#' @importFrom stats variable.names
#'
#' @export
#' 
variable.names.ser_normal <- function (object, ...)
  names(object$prior$prob)

# Compute the least-squares estimate and its variance given the data
# (X, y) and the residual variances (sigma).
compute_lse <- function (X, y, sigma) {
  xx   <- colSums(X^2)
  bhat <- drop(y %*% X)/xx
  shat <- sigma/xx
  names(bhat) <- colnames(X)
  names(shat) <- colnames(X)
  return(list(bhat = bhat,shat = shat))
}

# Compute the Bayes factor and posterior statistics for an EBNM model
# with a zero-centered normal prior. The inputs are the least-squares
# estimates (bhat), the variances of those estimates (s0), and the
# prior variance of the coefficients (s0). The outputs are the
# log-Bayes factor (lbf), the posterior mean estimates (b1), and the
# posterior variances (s1).
#
#' @importFrom stats dnorm
compute_posterior_ebnm_normal <- function (bhat, shat, s0) {
  lbf <- dnorm(bhat,0,sqrt(s0 + shat),log = TRUE) - 
         dnorm(bhat,0,sqrt(shat),log = TRUE)
  s1 <- 1/(1/s0 + 1/shat)
  b1 <- bhat*s1/shat
  return(list(lbf = lbf,b1 = b1,s1 = s1))
}

# Compute KL(p_post || p_ser), the Kullback-Leibler divergence between
# the SER posterior and the SER prior.
compute_kl_ser <- function (bhat, shat, b1, s1, lbf)
  -lbf + sum(b1*bhat/shat) - sum((s1 + b1^2)/shat)/2 
