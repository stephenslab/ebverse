#' @title Title goes here.
#'
#' @description Description goes here.
#' 
#' @export
#' 
prepare_susie_data <- function () {
  out <- list()
  class(out) <- c("susie_data","list")
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
prepare_susie_summary_data <- function () {
  out <- list()
  class(out) <- c("susie_data","list")
  return(out)
}
