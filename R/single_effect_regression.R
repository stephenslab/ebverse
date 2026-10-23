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
  fit <- list(prior = prior,posterior = as.numeric(NA),control = control)
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
       verbose = "none")

#' @rdname single_effect_regression
#'
#' @export
#'
fit_ser_model <- function (fit, dat) {
   UseMethod("fit_ser_model")
}

#' @rdname single_effect_regression
#'
#' @param dat Describe the "dat" input here.
#' 
#' @export
#' 
fit_ser_model.ser_normal <- function (fit, dat) {
  # TO DO: Implement estimation of the prior variance.
  w <- fit$prior$prob
  res <- compute_posterior_ebnm_normal(dat[,"bhat"],dat[,"s"],
                                       fit$prior$prior_variance)
  ans <- normalizelogweights(log(w) + res$lbf)
  pip <- ans$probs
  b1 <- pip * res$b1
  fit$posterior <- cbind(lbf = res$lbf,pip = pip,b1 = b1)
  fit$lbf <- ans$weight_sum
  return(fit)
}


