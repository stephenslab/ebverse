# This R file contains my very simple implementations of the single
# effect regression (SER) models.

# Compute the log-normalizing factor for the IG(a,b) distribution.
inv_gamma_factor <- function (a, b)
  a*log(b) - lgamma(a)

# Safely compute probabilities from log-weights.
normalizelogweights <- function (logw) {
  x <- max(logw)
  w <- exp(logw - x)
  return(w/sum(w))
}

# For each column of X, compute the least-squares estimate of the
# regression coefficient (bhat) and its variance (shat). Here, s is
# the residual variance.
least_squares_coefs <- function (X, y, s) {
  xx   <- colSums(X^2)
  xy   <- drop(crossprod(X,y))
  return(list(bhat = xy/xx,shat = s/xx))
}

# Compute the (log) Bayes factors given the least-squares estimates
# (bhat) and their corresponding variances (shat). Here, s0 is the prior
# variance of the regression coefficient.
compute_bayes_factors <- function (bhat, shat, s0)
  dnorm(bhat,0,sqrt(s0 + shat),log = TRUE) - 
  dnorm(bhat,0,sqrt(shat),log = TRUE)

# Compute the (log) Bayes factors and additional statistics under the
# normal-inverse-gamma (NIG) prior given the least-squares estimates
# (bhat) and their corresponding variances (shat). Additional outputs
# include the compute the posterior moments, mean (b1), second moment
# (b2) and variance (s1), and the posterior mode of the residual
# variance (rv). Here, s0 is the prior variance of the regression
# coefficient, n is the sample size, yy = sum(y^2), sxy = cor(x,y),
# and the residual variance has an IG(a0/2,b0/2) prior. I perform all
# these computations at once because it is just easier that way.
compute_stats_NIG <- function (n, xx, xy, yy, sxy, s0, a0, b0) {
  r0   <- s0/(s0 + 1/xx)
  rss  <- yy*(1 - r0*sxy^2)
  a1   <- a0 + n
  b1   <- b0 + rss  
  lbf  <- -(log(1 + s0*xx) + a1*log(b1/(b0 + yy)))/2
  # Note that the line above is the same as: 
  # lbf = (log(1 - r0) - a1*log(b1/(b0 + yy)))/2
  bhat <- xy/xx
  mu1  <- r0*bhat
  s1   <- b1/(a1 - 2)*r0/xx
  rv   <- (b1/2)/(a1/2 - 1)
  return(list(lbf = lbf,
              b1  = mu1,
              b2  = s1 + mu1^2,
              s1  = s1,
              rv  = rv))
}

# Compute the posterior moments: mean (b1), second moment (b2),
# variance (s1). The inputs are the least-squares estimates (bhat),
# their variances (shat), and the prior variance (s0).
compute_posterior_moments <- function (bhat, shat, s0) {
  s1 <- 1/(1/s0 + 1/shat)
  b1 <- s1*bhat/shat
  b2 <- b1^2 + s1
  return(list(b1 = b1,b2 = b2,s1 = s1))
}

# Compute the log-likelihood under the "null" model.
compute_null_loglik <- function (y, s)
  sum(dnorm(y,0,sqrt(s),log = TRUE))

# Compute the log-likelihood under the "null" model under the
# normal-inverse-gamma (NIG) prior.
compute_null_loglik_NIG <- function (n, yy, a0, b0) {
  return(-n*log(2*pi)/2 + 
         inv_gamma_factor(a0/2,b0/2) - 
         inv_gamma_factor((a0 + n)/2,(b0 + yy)/2))
}

# Compute the expected residual sum of squares. The inputs are the
# data (X, y); the posterior inclusion probabilities (pip), and the
# posterior second moments, the mean (b1) and the second moment (b2).
compute_erss <- function (X, y, pip, b1, b2) {
  xx  <- colSums(X^2)
  bb  <- pip * b1
  bb2 <- pip * b2
  xb  <- drop(X %*% bb)
  r   <- drop(y - xb)  
  return(sum(r^2) - sum(xb^2) + sum(xx * bb2))
}

# Compute the ELBO (which is also the marginal log-likelihood since
# the ELBO is exact). The inputs are: ll0, the log-likelihood under
# the "null" model; and lbf, the log-Bayes factors.
compute_elbo <- function (ll0, lbf) {
  x <- max(lbf)  
  bf <- exp(lbf - x)
  return(ll0 + x + log(mean(bf)))
}

# The EM update for the prior variance is the posterior second moment
# of the coefficient. The inputs are the posterior inclusion
# probabilities (pip) and the posterior second moments (b2).
update_prior_variance <- function (pip, b2)
  sum(pip * b2)

# Update the prior variance under the NIG prior by maximizing the
# marginal log-likelihood.
update_prior_variance_NIG_optim <- function (n, xx, xy, yy, sxy, 
                                             s0, a0, b0) {
  ll0 <- compute_null_loglik_NIG(n,yy,a0,b0)
  f <- function (s0) {
    lbf <- compute_stats_NIG(n,xx,xy,yy,sxy,s0,a0,b0)$lbf
    return(-compute_elbo(ll0,lbf))
  }
  out <- optim(s0,f,method = "Brent",lower = 0,upper = 10)
  return(out$par)
}

# This is an experimental implementation to check that my EM update is
# correct.
update_prior_variance_NIG_EM <- function (n, xx, xy, yy, sxy, pip, 
                                          s0, a0, b0) {
  p    <- length(pip)
  r0   <- s0/(s0 + 1/xx)
  rss  <- yy*(1 - r0*sxy^2)
  a1   <- a0 + n
  b1   <- b0 + rss  
  bhat <- xy/xx
  mu1  <- r0*bhat
  s1   <- r0/xx
  #
  # This was my previous (approximate) Monte Carlo implementation:
  #
  #   ns  <- 100000
  #   out <- 0
  #   for (i in 1:ns) {
  #     ss  <- 1/rgamma(p,a1/2,b1/2) # Same as: ss = rinvgamma(p,a1/2,b1/2)
  #     bb  <- rnorm(p,mu1/sqrt(ss),sqrt(s1))
  #     out <- out + sum(pip * bb^2)
  #   }
  #   return(out/ns)
  #
  u  <- gamma(1/2)/beta(a1/2,1/2)
  mb <- mu1 * sqrt(2/b1) * u
  vb <- s1 + mu1^2 * 2/b1 * (1/beta(a1/2,1) - u^2)
  return(sum(pip * (vb + mb^2)))
}

# My simple implementation of EM for the SER models. 
# For the "normal-inverse-gamma" (NIG) prior, in which the residual
# variance has an IG(a0/2,b0/2) prior. For this prior, note that
# prior mean = (b0/2)/(a0/2 - 1), 
# prior mode = (b0/2)/((a0/2) + 1)
simple_ser <- function (X, y, s0 = 1, s = 1, numiter = 10, 
                        prior = c("normal","NIG"), 
                        a0 = 0.01, b0 = 0.01, verbose = TRUE) {
  prior <- match.arg(prior)

  # Center X and y which is needed to account for an intercept.  
  X <- scale(X,center = TRUE,scale = FALSE)
  y <- drop(scale(y,center = TRUE,scale = FALSE))

  # Iterate the EM updates.
  n <- length(y)
  elbo <- rep(0,numiter)
  for (i in 1:numiter) {
    if (verbose)
      cat("(",i,") s=",s,", s0=",s0,"\n",sep = "")

    # Compute the Bayes factors and posterior inclusion probabilities
    # (PIPs).
    if (prior == "normal") {
      dat  <- least_squares_coefs(X,y,s)
      bhat <- dat$bhat
      shat <- dat$shat
      lbf  <- compute_bayes_factors(bhat,shat,s0)
    } else if (prior == "NIG") {
      dat <- list(xx  = colSums(X^2),
                  xy  = drop(crossprod(X,y)),
                  yy  = sum(y^2),
                  sxy = drop(cor(X,y)))
      lbf <- with(dat,compute_stats_NIG(n,xx,xy,yy,sxy,s0,a0,b0))$lbf
    }
    pip <- normalizelogweights(lbf)

    # Compute the ELBO.
    if (prior == "normal") {
      ll0 <- compute_null_loglik(y,s)
    } else if (prior == "NIG") {
      ll0 <- compute_null_loglik_NIG(n,dat$yy,a0,b0)
    }
    elbo[i] <- compute_elbo(ll0,lbf)

    # Compute the posterior means (b1), posterior variances (s1) and
    # posterior second moments (b2) of the coefficients.
    if (prior == "normal") {
      out <- compute_posterior_moments(bhat,shat,s0)
    } else if (prior == "NIG") {
      out <- with(dat,compute_stats_NIG(n,xx,xy,yy,sxy,s0,a0,b0))
      rv <- out$rv
    }
    b1 <- out$b1
    b2 <- out$b2
    s1 <- out$s1

    if (prior == "normal") {

      # The EM update for the residual variance is ERSS/n, where "ERSS" is
      # the expected residual sum of squares.
      erss <- compute_erss(X,y,pip,b1,b2)
      s <- erss/n
    } else if (prior == "NIG") {

      # Obtain the posterior mean of the residual variance parameter.
      s <- sum(pip * rv)
    }

    # The EM update for the prior variance is the posterior second
    # moment of the coefficient.
    if (prior == "normal") {
      s0 <- update_prior_variance(pip,b2)
    } else if (prior == "NIG") {

      # This update uses optimize().
      # s0 <- with(dat,
      #            update_prior_variance_NIG_optim(n,xx,xy,yy,sxy,s0,a0,b0))

      # This is the EM update.
      s0 <- with(dat,
                 update_prior_variance_NIG_EM(n,xx,xy,yy,sxy,pip,s0,a0,b0))
    }
  }

  return(list(s = s,s0 = s0,pip = pip,b1 = b1,b2 = b2,s1 = s1,elbo = elbo))
}

#' @rdname fit_susie_model
#'  
#' @title Title goes here.
#'
#' @param dat Describe "dat" here. (This should be a "prepared" data
#'   set.)
#' 
#' @param prior Describe "prior" here.
#'
#' @param L Describe "L" here.
#' 
#' @return Description of the output goes here.
#' 
#' @export
#' 
init_susie_model <- function (dat, prior, L = 10) {
  fit <- list()
  class(fit) <- c("susie","list")
  return(fit)
}

#' @title Title goes here.
#'
#' @description Description goes here.
#'
#' @param dat An object
#'
#' @seealso \code{\link{prepare_mlm_data}},
#'   \code{\link{prepare_mlm_summary_data}}
#'
#' @examples
#' # TO DO: Add example(s) here.
#' 
#' @export
#'
#' @keywords internal
#' 
fit_susie_model <- function (X, y, L = 10, s0 = rep(1,L), numiter = 10, 
                             a0 = 0.01, b0 = 0.01, verbose = TRUE,
                             control = list()) {

  # Initialize the posterior mean estimates of b, and initialize
  # storage for the PIPs and other quantities.
  n    <- nrow(X)
  p    <- ncol(X)
  B    <- matrix(0,p,L)
  PIP  <- matrix(0,p,L)
  elbo <- rep(0,L)

  # Initialize the expected residuals.
  r <- y

  # Iterate the SER updates.
  delta <- rep(0,numiter)
  for (i in 1:numiter) {
    if (verbose)
      cat(sprintf("iter %d: ",i))
    PIP0 <- PIP

    # Repeat for each single effect.
    for (k in 1:L) {
      b <- B[,k]
      r <- drop(r + X %*% b)
      out <- simple_ser(X,r,s0 = s0[k],numiter = 1,prior = "NIG", 
                        a0 = a0,b0 = b0,verbose = FALSE)
      b       <- with(out,pip * b1)
      B[,k]   <- b
      PIP[,k] <- out$pip
      s0[k]   <- out$s0
      elbo[k] <- out$elbo
      r <- drop(r - X %*% b)
    }

    # Record the largest change in the PIPs.
    delta[i] <- max(abs(PIP - PIP0))
    if (verbose)
      cat(sprintf("delta = %0.2g\n",delta[i]))
  }

  return(list(PIP = PIP,B = B,s0 = s0,elbo = elbo,delta = delta))
}

#' @rdname fit_susie_model
#' 
#' @export
#' 
fit_susie_control_default <- function ()
  list(verbose = "detailed")
