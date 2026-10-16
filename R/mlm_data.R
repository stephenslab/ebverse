#' @title Title goes here.
#'
#' @description Description goes here.
#' 
#' @export
#' 
prepare_mlm_data <- function () {
  out <- list()
  class(out) <- c("prepared_data_mlm","list")
  return(out)
}

#' @rdname prepare_mlm_data
#' 
#' @title Title goes here.
#'
#' @description Description goes here.
#' 
#' @export
#' 
prepare_mlm_summary_data <- function () {
  out <- list()
  class(out) <- c("prepared_data_mlm","list")
  return(out)
}
