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
#'   affect the model fitting. See \sQuote{Details} for more
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
#'   \item{lbf}{Approximate log-Bayes factor comparing the SuSiE model
#'     to the \dQuote{null} model (b = 0).}
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
#' @seealso \link{single_effect_regression}, \link{mlm_data}
#'
#' @details
#' TO DO: Explain that the special case of L = 1 reduces to
#' a single effect regression.
#' 
#' @examples
#' # TO DO: Add one or several examples.
#'
#' @importFrom Rcpp evalCpp
#' @importFrom Matrix sparseMatrix
#' @importFrom utils packageVersion
#' 
#' @useDynLib ebverse
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
    names(ser) <- paste0("L",1:L)
  control <- modifyList(susie_control_default(),control,keep.null = TRUE)
  J <- length(variable.names)
  ss <- data.frame(bhat = rep(as.numeric(NA),J),
                   shat = rep(as.numeric(NA),J),
                   xs   = rep(as.numeric(NA),J),
                   row.names = variable.names)
  fit <- list(elbo              = as.numeric(NA),
              lbf               = as.numeric(NA),
              ss                = ss,
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
#' @param dat An object of class \sQuote{mlm_data}, for example the
#'   output from a call to \code{\link{prepare_mlm_full_data}} or
#'   \code{\link{prepare_mlm_summary_data}}. See \link{mlm_data}
#'   for other ways to create \sQuote{mlm_data} objects.
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
  # TO DO NEXT: Implement various S3 methods: coef, plot, predict,
  #  summary, print logLik, residuals, variable.names
  # 
  # TO DO:
  # - Run R CMD check --as-cran.
  # - Write a function to simulate different X matrices---
  #   see https://github.com/stephenslab/dsc-linreg.
  # - Implement some functions for simulating data from susie models.
  # - Implement a "summary" S3 method for "ser_normal" objects.
  # - Implement a "summary" S3 method for "susie" objects.
  # - Implement "coef" S3 method.
  # - Implement "predict" S3 method.
  # - Implement a "plot" S3 method (idea: bhat vs. b1 scatterplot).
  # - Other S3 methods?
  # - Implement a new mlm_data object for PLINK "pgen" data sets.
  #   (Idea: specialized packages should be "Suggests".)
  # - Implement "global" estimation of parameters in the 
  #   SER models.
  #

  # Check and process the inputs.
  if (!inherits(fit,"susie"))
    stop("'fit' should be an object of class 'susie'",domain = NA)
  if (!inherits(dat,"mlm_data"))
    stop("'dat' should be an object of class 'mlm_data'",domain = NA)
  force(return_data)
  if (missing(dat)) {
    if (is.null(fit$dat))
      stop("the data were not provided",domain = NA)
    dat <- fit$dat
  }  
  if (any(attr(dat$ss,"residual_variance") != 1)) 
    stop("dat$ss does not give the correct summary statistics",domain = NA)
  if (any(variable.names(fit) != variable.names(dat)))
    stop("variables in data do not match variables declared in model")

  # Fit the SuSiE model and compute the posterior statistics.
  fit <- fit_susie_model_core(dat,fit,fit$control)

  # Prepare the output
  #
  # TO DO (e.g., add names to any outputs).
  #
  if (return_data)
    fit$dat <- dat
  else
    fit["dat"] <- NULL
  return(fit)
}

# This implements the core SuSiE model fitting algorithm.
fit_susie_model_core <- function (dat, fit, control) {
  progress <- as.numeric(NA)
  if (control$trace | control$verbose == "detailed") {
    progress <- as.data.frame(matrix(0,control$maxiter,6))
    names(progress) <- c("iter","elbo","resid.var","d.elbo","d.b","d.pip")
    progress[,"iter"] <- 1:control$maxiter
    if (control$verbose == "detailed") {
      cat("iter       ELBO   sigma  d(ELBO) delta(b)  d(pip) d(sigma)\n") 
      cat("---- ---------- ------- -------- -------- ------- --------\n")
    }
  }

  # Get the number of single effects.
  L <- length(fit$ser) 

  # Iterate the SER updates.
  for (iter in 1:control$maxiter) {

    # Save some quantities to monitor progress.
    b <- compute_post_mean_susie(fit$ser)
    current_state <- list(b    = b,
                          pip  = susie_get_pips(fit$ser),
                          elbo = fit$elbo,
                          residual_variance = fit$residual_variance)

    # Perform the SER updates.
    for (k in 1:L) {
      res <- susie_ibss_update(dat,fit$ser[[k]],L,b,fit$residual_variance)
      b <- res$b
      fit$ser[[k]] <- res$ser
    }

    # Compute the Bayes factor and ELBO.
    res      <- susie_compute_elbo(dat,fit)
    fit$lbf  <- res$lbf
    fit$elbo <- res$elbo
    erss     <- res$erss

    # Update the residual variance parameter.
    if (control$estimate_residual_variance)
      fit$residual_variance <- susie_update_residual_variance(dat,erss)

    # Update the progress data frame.
    if (is.data.frame(progress)) {
      pip <- susie_get_pips(fit$ser)
      progress[iter,"elbo"]        <- fit$elbo
      progress[iter,"resid.var"]   <- fit$residual_variance
      progress[iter,"d.elbo"]      <- fit$elbo - current_state$elbo
      progress[iter,"d.elbo"]      <- fit$elbo - current_state$elbo
      progress[iter,"d.b"]         <- max(abs(b - current_state$b))
      progress[iter,"d.pip"]       <- max(abs(pip - current_state$pip))
      progress[iter,"d.resid.var"] <- fit$residual_variance - 
                                      current_state$residual_variance
    }

    # Update the user about the algorithm's process.
    if (control$verbose == "detailed")
      cat(sprintf("%4d %+8.3e %07.4f %+8.1e %8.2e %7.1e %+7.1e\n",
                  progress[iter,"iter"],
                  progress[iter,"elbo"],
                  progress[iter,"resid.var"],
                  progress[iter,"d.elbo"],
                  progress[iter,"d.b"],
                  progress[iter,"d.pip"],
                  progress[iter,"d.resid.var"]))

    # Check for errors as well as the stopping criteria. Some checks
    # are only performed after the first iteration.
    if (fit$residual_variance <= 0)
      stop("prior variance is zero or negative",domain = NA)
    if (iter > 1) {
      if (!is.na(current_state$elbo) && 
          !is.na(control$elbotol) &&
          fit$elbo < current_state$elbo)
        warning("updates decreased the ELBO",domain = NA)
      if (!is.na(control$paramtol) &&
          max(abs(current_state$pip - pip)) < control$paramtol) {
        if (control$verbose != "none")
          message("parameter estimation tolerance is met; ",
                  "model fitting is stopped",
                  domain = NA)
        break
      }
      if (!is.na(current_state$elbo) && 
          !is.na(fit$elbo) &&
          !is.na(control$elbotol) &&
          fit$elbo - current_state$elbo < control$elbotol) {
        if (control$verbose != "none") 
          message("ELBO tolerance is met; model fitting is stopped",
                  domain = NA)
        break
      }
    }
  }

  # Prepare the output.
  fit$ss <- data.frame(bhat = dat$ss$bhat,
                       shat = dat$ss$shat,
                       xs   = dat$xs,
                       row.names = variable.names(dat))
  if (control$trace) {
    progress <- progress[1:iter,]
    if (is.data.frame(fit$progress)) {
      progress$iter <- progress$iter + max(fit$progress$iter)
      fit$progress <- rbind(fit$progress,progress)
    }
    else 
      fit$progress <- progress
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

#' @importFrom stats variable.names
#'
#' @export
#' 
variable.names.susie <- function (object, ...)
  rownames(object$ss)

# Compute the posterior inclusion probabilities (PIPs) given the
# single effect regression PIPs.
susie_get_pips <- function (ser_list) {
  varnames <- variable.names(ser_list[[1]])
  J <- length(varnames)
  L <- length(ser_list)
  pip <- rep(1,J)
  for (k in 1:L)
    pip <- pip * (1 - ser_list[[k]]$posterior$pip)
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
                              shat = setNames(sigma*d,varnames),
                              residual_variance = 1)
  ser <- fit_ser_model(ser,dat_ser)
  bk <- ser$posterior$mean
  sk <- ser$posterior$var
  if (L > 1)
    b <- b + bk

  return(list(ser = ser,b = b))
}

# Compute the log-likelihood under the "null" (b = 0) model.
compute_null_loglik <- function (dat, sigma) {
  n <- dat$n
  yy <- (n-1)*dat$var_y
  return(-(yy/sigma + n*log(2*pi*sigma))/2)
}

# Compute the "expected sum of squares" (ERSS) for the SuSiE model,
# omitting the crossprod(y) term.
susie_compute_erss_without_yy <- function (dat, ser_list) {
  L <- length(ser_list)
  b <- compute_post_mean_susie(ser_list)
  d <- dat$ss$shat
  xy <- dat$ss$bhat/d
  erss <- sum(b*compute_XtXb(dat,b)) - 2*sum(b*xy)
  for (k in 1:L) {
    bk   <- ser_list[[k]]$posterior$mean
    sk   <- ser_list[[k]]$posterior$var
    erss <- erss - sum(bk*compute_XtXb(dat,bk)) + sum((sk + bk^2)/d)
  }
  return(erss)
}

# Return (1) ERSS, (2) approximate log-Bayes factor, (3) ELBO. See
# the susie_compute_erss_without_yy function for what we mean by the
# "ERSS".
susie_compute_elbo <- function (dat, fit) {
  L     <- length(fit$ser)
  sigma <- fit$residual_variance
  llik0 <- compute_null_loglik(dat,sigma)
  erss  <- susie_compute_erss_without_yy(dat,fit$ser)
  if (L == 1) {
    KL  <- fit$ser[[1]]$KL
    lbf <- fit$ser[[1]]$lbf
  } else {
    KL  <- sapply(fit$ser,"[[","KL")
    lbf <- -erss/(2*sigma) - sum(KL)
  }
  return(list(erss = erss,lbf = lbf,elbo = llik0 + lbf))
}

# Output the EM update for the residual variance. See the
# susie_compute_erss_without_yy function for an explanation of the
# "erss" input.
susie_update_residual_variance <- function (dat, erss) {
  n <- dat$n
  yy <- (n-1)*dat$var_y
  return((yy + erss)/n)
}
