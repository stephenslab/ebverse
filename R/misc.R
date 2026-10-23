# Given a vector of log-weights (lw), compute the following in a
# numerically stafe way: (1) weight_sum = sum(w), and (2) probs =
# w/sum(w), where w = exp(lw).
normalizelogweights <- function (lw) {
  a <- max(lw)
  w <- exp(lw - a)
  return(list(weight_sum = log(sum(w)) + a,probs = w/sum(w)))
}
