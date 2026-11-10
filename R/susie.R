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
#' @param variable.names The names of the variables in the SuSiE
#'   model.
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
  function (variable.names, 
            L = 10, 
            ser = init_ser_model(variable.names,family = "normal",
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
    stop("L should be at least 1",domain = NA)
  if (inherits(ser,"ser")) {
    ser <- rep(list(ser),L)
  } else if (is.list(ser)) {
    if (length(ser) != L)
      stop("'ser' does not match 'L'",domain = NA)
  } else
    stop("'ser' should be an object of class 'ser', or a list of ser ",
         "objects",domain = NA)
  if (is.null(names(ser)))
    names(ser) <- paste0("single_effect_",1:L)
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
  # - Write a function to simulate different X matrices---
  #   see https://github.com/stephenslab/dsc-linreg.
  # - Implement some functions for simulating data from susie models.
  # - Test fit_susie_model with L = 1.
  # - Test *refitting* of a susie model.
  # - Implement the "susie_data" functions.
  # - Implement estimation of the residual variance.
  # - Implement the ELBO computations.
  # - Implement a "summary" S3 method.
  # - Implement "global" estimation of parameters in the 
  #   SER models.
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
  if (control$trace | control$verbose == "detailed") {
    progress <- as.data.frame(matrix(0,control$maxiter,5))
    names(progress) <- c("iter","elbo","resid.var","d.b","d.pip")
    progress[,"iter"] <- 1:control$maxiter
    if (control$verbose == "detailed") {
      cat("iter      ELBO sigma delta(b)  d(pip)\n") 
      cat("---- --------- ----- -------- -------\n")
    }
  }

  # Get the number of single effects.
  L <- length(fit$ser) 

  # Iterate the SER updates.
  for (iter in 1:control$maxiter) {

    # Compute the posterior mean of b and the PIPs.
    b   <- compute_post_mean_susie(fit$ser)
    pip <- susie_get_pips(fit)

    # Save some quantities to monitor progress.
    current_state <- list(b = b,
                          pip = pip,
                          residual_variance = fit$residual_variance)

    # Repeat for each single effect.
    for (k in 1:L) {

      # Perform an SER update.
      res <- susie_ibss_update(dat,fit$ser[[k]],L,b,fit$residual_variance)
      b <- res$b
      fit$ser[[k]] <- res$ser
    }

    # Update the residual variance parameter.
    if (control$estimate_residual_variance) {
      fit$residual_variance <- susie_update_residual_variance(dat,fit$ser)
    }

    #
    # TO DO
    #
    # The EM update for the residual variance is ERSS/n, where "ERSS" is
    # the expected residual sum of squares.
    #  erss <- compute_erss(X,y,pip,b1,b2)
    #  s <- erss/n
    #

    # Update the progress data frame.
    if (control$trace | control$verbose == "detailed") {
      pip <- susie_get_pips(fit)
      progress[iter,"elbo"]        <- 0
      progress[iter,"resid.var"]   <- fit$residual_variance
      progress[iter,"d.b"]         <- max(abs(b - current_state$b))
      progress[iter,"d.pip"]       <- max(abs(pip - current_state$pip))
      progress[iter,"d.resid.var"] <- fit$residual_variance - 
                                      current_state$residual_variance
    }

    # Update the user about the algorithm's process.
    if (control$verbose == "detailed")
      cat(sprintf("%4d %8.3e %5.2f %8.2e %7.1e\n",
                  progress[iter,"iter"],
                  progress[iter,"elbo"],
                  progress[iter,"resid.var"],
                  progress[iter,"d.b"],
                  progress[iter,"d.pip"]))

      # Check for errors as well as the stopping criteria.
    if (fit$residual_variance <= 0)
      stop("prior variance is zero or negative",domain = NA)
    # TO DO: Check for decreases in the ELBO.
    if (!is.na(control$paramtol) &&
        max(abs(current_state$pip - pip)) < control$paramtol) {
      if (control$verbose != "none")
        message("parameter estimation tolerance is met; ",
                "model fitting is stopped",
                domain = NA)
      break
    }
    # TO DO: Check whether the ELBO increases by a smaller amount than
    # the tolerance.
  }

  # Prepare the output.
  if (control$trace) {
    progress  <- progress[1:iter,]
    if (is.data.frame(fit$progress)) {
      progress$iter <- progress$iter + max(fit$progress$iter)
      fit$progress  <- rbind(fit$progress,progress)
    }
    else 
      fit$progress <- progress
    rownames(fit$progress) <- NULL
  } else
    fit$progress <- as.numeric(NA)
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
       trace    = TRUE,
       verbose  = "detailed")

# TO DO: Explain here what this function does, and how to use it.
susie_get_pips <- function (fit) {
  varnames <- variable.names(fit$ser[[1]])
  J <- length(varnames)
  L <- length(fit$ser)
  pip <- rep(1,J)
  for (k in 1:L)
    pip <- pip * (1 - fit$ser[[k]]$posterior$pip)
  pip <- 1 - pip
  names(pip) <- varnames
  return(pip)
}

# Compute the posterior mean of "b" in a SuSiE model using the
# posterior means of the coefficients in the individual single effect
# regressions. Input "ser" is a list of single effect regressions.
compute_post_mean_susie <- function (ser_list) {
  L <- length(ser_list)
  varnames <- variable.names(ser_list[[1]])
  J <- length(varnames)
  b <- rep(0,J)
  for (k in 1:L)
    b <- b + ser_list[[k]]$posterior$mean
  names(b) <- varnames
  return(b)
}

# This implements the core update of the IBSS algorithm for susie,
# which involves updating the (approximate) posterior distribution for
# one SER while the other SERs remain fixed. Here, "ser" the SER that
# will be updated.
#
#' @importFrom stats setNames
susie_ibss_update <- function (dat, ser, L, b, sigma) {
  bhat     <- dat$ss$bhat
  d        <- dat$ss$shat
  varnames <- variable.names(ser)
  if (L > 1) {
    b    <- b - ser$posterior$mean
    r    <- compute_XtXb(dat,b)
    bhat <- bhat - d*r
  } 
  dat_ser <- prepare_ser_data(bhat = setNames(bhat,varnames),
                              shat = setNames(sigma*d,varnames))
  ser <- fit_ser_model(ser,dat_ser)
  if (L > 1)
    b <- b + ser$posterior$mean
  return(list(ser = ser,b = b))
}

# TO DO: Explain here what this function does and how to use it.
susie_update_residual_variance <- function (dat, ser) {
  B <- s$alpha * s$mu
  XB2 <- sum((B %*% XtX) * B)
  betabar <- colSums(B)
  d <- attr(XtX,"d")
  postb2 <- s$alpha * s$mu2 # Posterior second moment.
  ER2_ss <- yty - 2*sum(betabar * Xty) + 
  sum(betabar * (XtX %*% betabar)) -
  XB2 + sum(d * t(postb2)))
  fit$residual_variance <- ER2_ss/n
}
