#' @title Sample Names
#'
#' @description Simple utility for extacting the names of the samples
#'   from a data object.
#'
#' @param dat A data object.
#'
#' @param \dots Other arguments.
#'
#' @return A character vector containing the names of the samples.
#' 
#' @export
#'
sample.names <- function (dat, ...) {
  UseMethod("sample.names")
}

