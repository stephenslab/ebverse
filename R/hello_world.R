#' @useDynLib ebverse
#'
#' @export
#'
hello_world <- function (x) {
  cat("Hello world!\n")
  return(rcpp_hello_world());
}
