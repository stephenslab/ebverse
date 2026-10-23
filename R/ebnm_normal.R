# TO DO: Explain here what this function is for, and what are the
# inputs and outputs.
#
#' @importFrom stats dnorm
compute_posterior_ebnm_normal <- function (bhat, s, s0) {
  lbf <- dnorm(bhat,0,sqrt(s0 + s),log = TRUE) - 
         dnorm(bhat,0,sqrt(s),log = TRUE)
  s1 <- 1/(1/s0 + 1/s)
  b1 <- bhat*s1/s
  return(list(lbf = lbf,b1 = b1,s1 = s1))
}
