#' @name mlm_data
#'
#' @title Data Objects for Multiple Linear Regression Models
#'
#' @description Functions to prepare various types of data for 
#'
#' @param X An n x J matrix containing observations about J candidate
#'   explanatory variables.
#'
#' @param y A vector of regression outcomes.
#'
#' @param standardize If \code{standardize = TRUE}, columns of X are
#'   standardized. If variables are in the same units already, you
#'   might not wish to standardize. Coefficients are always returned
#'   on the original scale.
#'
#' @return An \sQuote{mlm_data} object.
#'
#' @details Add details here.
#'
#'   \subsection{Design of MLM data classes}{
#'     TO DO. 
#'   }
#' 
#' @export
#'
prepare_mlm_full_data <- function (X, y, standardize = FALSE) {
  #
  # TO DO:
  #
  # - Define "mlm_data" class in the documentation.
  # - Implement intercept = TRUE.
  # - Add intercept estimate to data object.
  #

  # Check and process the inputs.  
  if (!(is.numeric(y) & 
        length(y) > 1))
    stop("y should be a numeric vector containing at least 2 data points",
         domain = NA)
  if (!(is.numeric(X) &
        is.matrix(X) &
        nrow(X) == length(y) &
        nrow(X) > 1 &
        ncol(X) > 1))
    stop("X should be a numeric matrix with at least 2 rows and 2 ",
         "columns, and nrow(X) should be the same as length(y)",
         domain = NA)

  # Check the variable and sample names.
  if (is.null(names(y)) | 
      is.null(rownames(X)) |
      is.null(colnames(X)))
    stop("Elements of y, and rows and columns of X, need to be named",
         domain = NA)
  if (any(names(y) != rownames(X)))
    stop("Sample names in X and y are not consistent",domain = NA)

  # For now, automatically center X and y. 
  # (TO DO: make this an option).
  # 
  # Standardize the columns of X, if requested.
  if (standardize) {
    X <- scale(X,center = TRUE,scale = TRUE)
    xs <- attr(X,"scaled:scale")
  } else {
    X <- scale(X,center = TRUE,scale = FALSE)
    xs <- rep(1,J)
  }

  # Prepare the output.
  n <- nrow(X)
  J <- ncol(X)
  dat <- list(n     = n,
              nvar  = J,
              X     = X,
              y     = y,
              var_y = sum(y^2)/(n-1),
              ss    = prepare_ser_data(X = X,y = y),
              xs    = xs,
              standardize = standardize)
  class(dat) <- c("mlm_full_data","mlm_data","list")
  return(dat)
}

#' @rdname mlm_data
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
prepare_mlm_summary_data <- function (z, R, bhat, shat, n, var_y) {
  #  
  # NOTES: 
  # - Be strict about naming things.
  # - I need to do a better job checking the inputs.
  # - Handle situation when n and/or var_y are missing.
  # - Give messages/warnings to tell user what is going on.
  # - See susie_rss.R.
  # - Assume n is provided, but not necessarily var_y.
  #
  # TO DO:
  # - Check that susie with z-scores produces the same result as
  #   susie with full data when X is standardized. (I don't think 
  #   Y needs to be standardized so long as var_y is provided?)
  #
  dat <- list()
  class(dat) <- c("mlm_summary_data","mlm_data","list")
  return(dat)
}

#' @title Sample Names
#'
#' @description Simple utility for extacting the names of the samples
#'   from a data object.
#'
#' @param dat A data object.
#'
#' @return A character vector containing the names of the samples.
#' 
#' @export
#'
sample.names <- function (dat, ...) {
  UseMethod("sample.names")
}

#'
#' @export 
#' 
sample.names.mlm_full_data <- function (object, ...)
  names(object$y)

#'
#' @importFrom stats variable.names
#' 
#' @export  
#'
variable.names.mlm_full_data <- function (object, ...)
  colnames(object$X)

#' 
#' @export
#' 
summary.mlm_full_data <- function (object, ...) {
  out <- list(n = object$n,
              nvar = object$nvar,
              standardize = object$standardize)
  class(out) <- c("summary_mlm_full_data","list")
  return(out)
}

#'
#' @export
#' 
print.summary_mlm_full_data <- function (x, ...) {
  cat(sprintf(paste("Sample-level or individual-level data for %d variables",
                    "in %d samples.\n"),x$n,x$nvar))
  cat(sprintf("Variables %s standardized.\n",
              ifelse(x$standardize,"are","are not")))
  return(invisible(x))
}

#'
#' @export
#' 
print.mlm_full_data <- function (x, ...) {
  print(summary(x))
  return(invisible(x))
}

#' @rdname mlm_data_internal
#'
#' @title Add Title Here.
#'
#' @description Add description here.
#'
#' @param dat Describe input argument "data" here.
#'
#' @param b Describe input argument "b" here.
#' 
#' @return Describe the return value here.
#' 
#' @seealso \link{mlm_data}
#' 
#' @export
#'
#' @keywords internal
#' 
compute_XtXb <- function (dat, b) {
  UseMethod("compute_XtXb")
}

#' @importFrom stats setNames
#' 
#' @export
#' 
compute_XtXb.mlm_full_data <- function (dat, b)
  setNames(drop(t(dat$X %*% b) %*% dat$X),colnames(dat$X))

