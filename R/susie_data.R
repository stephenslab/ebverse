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
  # - Be strict about naming things.
  # - I need to do a better job checking the inputs.
  # - Give messages/warnings to tell user what is going on.
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
#' @param var_y Describe input argument "yty" here.
#' 
#' @export
#' 
prepare_susie_summary_data <- function (z, R, bhat, shat, n, var_y) {
  #  
  # NOTES: 
  # - Be strict about naming things.
  # - I need to do a better job checking the inputs.
  # - Give messages/warnings to tell user what is going on.
  #
  dat <- list()
  class(dat) <- c("susie_summary_data","susie_data","list")
  return(dat)
}

#
# TO DO: Define the S3 methods for susie_data objects.
#
