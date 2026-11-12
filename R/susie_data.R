#' @name susie_data
#'
#' @title Title goes here.
#'
#' @description Description goes here.
#'
#' @param X Describe input argument "X" here.
#'
#' @param y Describe input argument "y" here.
#'
#' @seealso \link{susie}
#'
#' @export
#'
prepare_susie_full_data <- function (X, y) {
  #
  # NOTES: 
  # - Be strict about naming things!
  # - I need to do a better job checking the inputs.
  # - Give messages/warnings to tell user what is going on.
  #
  # TO DO:
  # - Implement variable.names S3 method.
  # - Implement sample.names S3 method.
  # - Get this working for both dense and sparse matrices.
  #   For sparse X, the output will be a "susie_full_data_sparse" 
  #   object.
  #
  dat <- list()
  class(dat) <- c("susie_full_data","susie_data","list")
  return(dat)
}

#' @rdname susie_data
#' 
#' @title Title goes here.
#'
#' @description Description goes here.
#'
#' @param z Describe input argument "XtX" here.
#'
#' @param R Describe input argument "R" here.
#'
#' @param bhat Describe input argument "bhat" here.
#'
#' @param shat Describe input argument "shat" here.
#' 
#' @param n Describe input argument "n" here.
#'
#' @param var_y Describe input argument "var_y" here.
#' 
#' @export
#' 
prepare_susie_summary_data <- function (z, R, bhat, shat, n, var_y) {
  #  
  # NOTES: 
  # - Be strict about naming things.
  # - I need to do a better job checking the inputs.
  # - Give messages/warnings to tell user what is going on (?).
  #
  dat <- list()
  class(dat) <- c("susie_summary_data","susie_data","list")
  return(dat)
}

#' @title Add Title Here.
#'
#' @description Add description here.
#'
#' @param dat Describe input argument "data" here.
#'
#' @return Describe the return value here.
#' 
#' @seealso \link{susie_data}
#' 
#' @export
#'
#' @keywords internal
#' 
compute_XtXb <- function (dat, b) {
  UseMethod("compute_XtXb")
}

#' 
#' @export
#' 
# TO DO: Update this function to accept either a matrix or a vector.
compute_XtXb.susie_full_data <- function (dat, b) {
  Y <- t(dat$X %*% b) %*% dat$X
  colnames(Y) <- colnames(dat$X)
  return(drop(Y))
}
