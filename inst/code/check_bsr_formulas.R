# Some code to check my equations for the Bayesian simple regression
# model.
set.seed(1)

# Check that: (i) likelihood under model y ~ N(xb,sigma^2*I) is the
# same as the likelihood under the EBNM model, bhat ~ N(b,s^2), up to
# a constant of proportionality; and (ii) the two expressions for the
# likelihood ratios are the same.
n <- 20
s0 <- 2
sigma <- 0.3
b <- 0.5
x <- rnorm(n)
x <- x - mean(x)
e <- rnorm(n,sd = sqrt(sigma))
y <- x*b + e
y <- y - mean(y)
b_true <- b
b <- seq(-4,4,length.out = 100)
ns <- length(b)
xx <- sum(x^2)
xy <- sum(x*y)
bhat <- xy/xx
s <- sqrt(sigma/xx)
z <- bhat/s
ll <- rep(0,ns)
lbf <- rep(0,ns)
ll_nm <- rep(0,ns)
lbf_nm <- rep(0,ns)
for (i in 1:ns) {
  ll[i] <- sum(dnorm(y,x*b[i],sqrt(sigma),log = TRUE))
  lbf[i] <- ll[i] - sum(dnorm(y,0,sqrt(sigma),log = TRUE))
  ll_nm[i] <- dnorm(bhat,b[i],s,log = TRUE)
  lbf_nm[i] <- ll_nm[i] - dnorm(bhat,0,s,log = TRUE)
}
d <- max(ll - ll_nm)
par(mfrow = c(1,3))
plot(b,ll,type = "l",col = "magenta",xlab = "b",ylab = "loglik")
lines(b,ll_nm + d,col = "black",lty = "dotted")
plot(b,lbf,type = "l",col = "magenta",xlab = "b",ylab = "logBF")
lines(b,lbf_nm,col = "black",lty = "dotted")
print(range(lbf - lbf_nm))

# Check my expression for the Bayes factor for the EBNM model with a
# zero-centered normal prior and the expression (A.3) in the susie
# paper (Wang et al, 2020).
n   <- 1000
bhat <- rnorm(n)
s0   <- abs(rnorm(n))
s    <- abs(rnorm(n))
z    <- bhat/sqrt(s)
lbf <- dnorm(bhat,0,sqrt(s + s0),log = TRUE) - 
              dnorm(bhat,0,sqrt(s),log = TRUE)
lbf_wang <- log(s/(s0 + s))/2 + z^2*s0/(s0 + s)/2
plot(lbf,lbf_wang,pch = 20)
abline(a = 0,b = 1,col = "magenta",lty = "dotted")
print(range(lbf - lbf_wang))
