# Compute the log-likelihood under the "null" model.
compute_null_loglik <- function (y, s)
  sum(dnorm(y,0,sqrt(s),log = TRUE))

# Compute the expected residual sum of squares. The inputs are the
# data (X, y); the posterior inclusion probabilities (pip), and the
# posterior second moments, the mean (b1) and the second moment (b2).
compute_erss <- function (X, y, pip, b1, b2) {
  xx  <- colSums(X^2)
  bb  <- pip * b1
  bb2 <- pip * b2
  xb  <- drop(X %*% bb)
  r   <- drop(y - xb)  
  return(sum(r^2) - sum(xb^2) + sum(xx * bb2))
}

# Compute the ELBO (which is also the marginal log-likelihood since
# the ELBO is exact). The inputs are: ll0, the log-likelihood under
# the "null" model; and lbf, the log-Bayes factors.
compute_elbo <- function (ll0, lbf) {
  x <- max(lbf)  
  bf <- exp(lbf - x)
  return(ll0 + x + log(mean(bf)))
}

#' @name susie
#'  
#' @title Fitting Sum of Single Effect (SuSiE) Models
#'
#' @description Methods for fitting \dQuote{Sum of Single Effect}
#'   (SuSiE) models and computing posterior statistics from SuSiE
#'   models. The \sQuote{fit_susie_model} method is the main method
#'   for fitting SuSiE models and computing posterior statistics.
#'
#' @param nvar The number of variables in the SuSiE model.
#' 
#' @param L An integer, 1 or greater, specifying the number of single
#'   effects.
#'
#' @param ser An object of class \sQuote{ser} or a list of
#'   \sQuote{ser} objects specifying the SER family for each of the
#'   L single effect regressions. If \code{ser} is a list, it should
#'   have L elements.
#'
#' @param residual_variance Initial estimate of the residual variance
#'   parameter.
#'
#' @param control An object, possibly a list, specifying settings that
#'   affect the model fitting. See \sQuote{details} for further
#'   information.
#'  
#' @return \code{init_susie_model} and \code{fit_susie_model} return
#'   objects of class \sQuote{susie} storing key posterior statistics
#'   for all variables and the state of the SuSiE model fit. It is a
#'   list with the following elements:
#'
#'   \item{ser}{The L single effect regression models which define the
#'     SuSiE posterior distribution.}
#' 
#'   \item{residual_variance}{The residual variance parameter.}
#'
#'   \item{elbo}{The value of the ELBO (an approximate marginal
#'     likelihood) at the current model fit.}
#' 
#'   \item{progress}{A data frame containing detailed information about
#'     progress on the model fitting.}
#'
#'   \item{control}{A list containing the settings that affect how the
#'     model fitting is performed.}
#' 
#'   \item{version}{The version of the ebverse R package that was
#'     used to generate this object}
#' 
#'   \item{dat}{When \code{return_data = TRUE}, this contains the 
#'     data that were used to fit the model.}
#'
#' @seealso \link{susie_data}, \link{single_effect_regression}
#'
#' @examples
#' # TO DO: Add one or several examples.
#'
#' @importFrom utils packageVersion
#' 
#' @export
#' 
init_susie_model <- 
  function (nvar, 
            L = 10, 
            ser = init_ser_model(nvar,family = "normal",
                                 control = list(maxiter = 1,trace = FALSE)),
            residual_variance = 1, 
            control = list()) {
  #
  # NOTES: 
  # - This is a function that is intended to be part of the main
  #   user interface and so there should be more careful checking
  #   of the inputs.
  # 
  if (L < 1)
    stop("L should be at least 1")
  if (inherits(ser,"ser")) {
    ser <- rep(list(ser),L)
  } else if (is.list(ser)) {
    if (length(ser) != L)
      stop("'ser' does not match 'L'")
  } else
    stop("'ser' should be an object of class 'ser', or a list of ser ",
         "objects",domain = NA)
  control <- modifyList(susie_control_default(),control,keep.null = TRUE)
  fit <- list(elbo              = as.numeric(NA),
              ser               = ser,
              residual_variance = residual_variance,
              progress          = as.numeric(NA),
              control           = control)
  fit$version <- packageVersion("ebverse")
  class(fit) <- c("susie","list")
  return(fit)
}

#' @rdname susie
#'
#' @param fit An object of class \sQuote{susie} that captures the
#'   state of a SuSiE model fit, typically the output from running
#'   an \sQuote{init_susie_model} function, or an output from a
#'   previous \code{fit_susie_model} call.
#'
#' @param dat An object of class \sQuote{susie_data}, for example the
#'   output from a call to \code{\link{prepare_susie_full_data}} or
#'   \code{\link{prepare_susie_summary_data}}. See \link{susie_data}
#'   for other ways to create \sQuote{susie_data} objects.
#' 
#' @param return_data If \code{return_data = TRUE}, the input data are
#'   included in the output.
#' 
#' @examples
#' # TO DO: Add example(s) here.
#' 
#' @export
#'
fit_susie_model <- function (fit, dat, return_data = missing(dat)) {
  #
  # NOTES: 
  # - This is a function that is intended to be part of the main
  #   user interface and so there should be more careful checking
  #   of the inputs.
  # 
  # TO DO:
  # - Test fit_susie_model with L = 1.
  # - Test refitting of a susie model.
  # - Add a "trace" control parameter (default: TRUE).
  # - Add a "progress" output.
  # - Implement the "susie_data" functions.
  # - Implement estimation of the residual variance.
  # - Implement the ELBO computations.
  # - Implement a "summary" S3 method.
  # - Add names to outputs wherever possible.
  #

  # Check and process the inputs.
  if (!inherits(fit,"susie"))
    stop("'fit' should be an object of class 'susie'",domain = NA)
  if (!inherits(dat,"susie_data"))
    stop("'dat' should be an object of class 'susie_data'",domain = NA)
  force(return_data)
  if (missing(dat)) {
    if (is.null(fit$dat))
      stop("the data were not provided",domain = NA)
    dat <- fit$dat
  }  
  if (any(attr(dat$ss,"residual_variance") != 1)) 
    stop("dat$ss does not give the correct summary statistics",domain = NA)

  # Fit the SuSiE model and compute the posterior statistics.
  fit <- fit_susie_model_core(dat,fit,fit$control)

  # Prepare the output
  #
  # TO DO
  #
  if (return_data)
    fit$dat <- dat
  else
    fit["dat"] <- NULL
  return(fit)
}

# This implements the core SuSiE model fitting algorithm.
fit_susie_model_core <- function (dat, fit, control) {

  # Get the number of single effects.
  L <- length(fit$ser) 

  # Iterate the SER updates.
  for (i in 1:control$maxiter) {

    # Compute the posterior mean of b.
    #
    # TO DO: Hide this code inside a separate (internal) function
    # called "compute_post_mean_susie".
    # 
    if (L > 1) {
      b <- rep(0,dat$nvar)
      for (k in 1:L)
        b <- b + fit$ser[[k]]$posterior$mean
    }
    current_state <- list(b = b)

    # Repeat for each single effect.
    for (k in 1:L) {

      # Perform an SER update.
      # 
      # TO DO: Hide parts of this code in a separate (internal)
      # function called "ser_update".
      # 
      fit_ser <- fit$ser[[k]]
      bhat    <- dat$ss$bhat
      d       <- dat$ss$shat
      sigma   <- fit$residual_variance
      if (L > 1) {
        bk   <- fit_ser$posterior$mean
        b    <- b - bk
        r    <- compute_XtXb(dat,b)
        bhat <- bhat - d*r
      } 
      dat_ser <- prepare_ser_data(bhat = bhat,shat = sigma*d)
      fit_ser <- fit_ser_model(fit_ser,dat_ser)
      if (L > 1) {
        bk <- fit_ser$posterior$mean
        b  <- b + bk
      }
      fit$ser[[k]] <- fit_ser
    }

    #
    # TO DO
    #
    # The EM update for the residual variance is ERSS/n, where "ERSS" is
    # the expected residual sum of squares.
    #  erss <- compute_erss(X,y,pip,b1,b2)
    #  s <- erss/n
    #

    #
    # TEMPORARY
    #
    cat(max(abs(current_state$b - b)),"\n")
  }

  return(fit)
}
 
#' @rdname susie
#' 
#' @export
#' 
susie_control_default <- function() 
  list(estimate_residual_variance = TRUE,
       paramtol = 1e-6,
       elbotol  = as.numeric(NA),
       maxiter  = 100,
       version  = "R",
       verbose  = "detailed")
