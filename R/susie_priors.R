#' @rdname susie_priors
#'
#' @title Title goes here.
#'
#' @description Description goes here.
#' 
#' @export
#' 
init_susie_prior_normal <- function () {
  out <- list()
  class(out) <- c("susie_prior_normal","susie_prior","eb_prior","list")
  return(out)
}

#' @rdname susie_priors
#'
#' @title Title goes here.
#'
#' @description Description goes here.
#' 
#' @export
#' 
init_susie_prior_laplace <- function () {

}
