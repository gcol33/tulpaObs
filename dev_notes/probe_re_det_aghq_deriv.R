# Finite-difference verification of the detection-arm AGHQ per-site marginal
# derivatives (R/re_aghq.R, arm == "det" branch of site_deriv). A bug in d1/d2
# would silently bias the per-group mode / curvature, so check them against
# numerical derivatives of logL before trusting the end-to-end recovery.
#
#   p = sigmoid(eta), psi fixed
#   detected (d > 0):  logL = log psi + d log p + (n-d) log(1-p)
#   non-detected (d=0): logL = log(psi (1-p)^n + (1 - psi))

cl <- function(e) pmin(pmax(e, -30), 30)

# Analytic d1/d2 transcribed verbatim from the arm == "det" site_deriv.
analytic <- function(eta, n, d, psi, detected) {
  p <- plogis(eta)
  if (detected) {
    logL <- log(psi) + d * log(p) + (n - d) * log1p(-p)
    d1 <- d - n * p
    d2 <- -n * p * (1 - p)
  } else {
    q <- exp(n * log1p(-p)); L <- psi * q + (1 - psi)
    A  <- -psi * n * p * q
    dA <- -psi * n * p * q * ((1 - p) - n * p)
    logL <- log(L)
    d1 <- A / L
    d2 <- (dA * L - A^2) / L^2
  }
  list(logL = logL, d1 = d1, d2 = d2)
}

logL_only <- function(eta, n, d, psi, detected) {
  p <- plogis(eta)
  if (detected) log(psi) + d * log(p) + (n - d) * log1p(-p)
  else log(psi * exp(n * log1p(-p)) + (1 - psi))
}

fd <- function(f, x, h = 1e-5) (f(x + h) - f(x - h)) / (2 * h)
fd2 <- function(f, x, h = 1e-4) (f(x + h) - 2 * f(x) + f(x - h)) / h^2

set.seed(7)
maxerr1 <- 0; maxerr2 <- 0
for (rep in 1:2000) {
  n   <- sample(1:6, 1)
  d   <- if (runif(1) < 0.5) 0L else sample(1:n, 1)
  psi <- runif(1, 0.05, 0.95)
  eta <- runif(1, -4, 4)
  detected <- d > 0
  a  <- analytic(eta, n, d, psi, detected)
  f  <- function(e) logL_only(e, n, d, psi, detected)
  e1 <- abs(a$d1 - fd(f, eta))
  e2 <- abs(a$d2 - fd2(f, eta))
  maxerr1 <- max(maxerr1, e1)
  maxerr2 <- max(maxerr2, e2)
}
cat(sprintf("max |d1 - FD|  = %.3e\n", maxerr1))
cat(sprintf("max |d2 - FD2| = %.3e\n", maxerr2))
cat(if (maxerr1 < 1e-4 && maxerr2 < 1e-2) "DERIVATIVES OK\n" else "DERIVATIVE MISMATCH\n")
