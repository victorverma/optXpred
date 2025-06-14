# ChatGPT R Translation of rstab.m MATLAB package

# ---- d2 function ----
d2 <- function(z) {
  p <- c(840.066852536483239, 20.0011141589964569)
  q <- c(1680.13370507296648, 180.0133704073900023, 0.1)
  y <- numeric(length(z))
  
  I1 <- abs(z) > 0.1
  I2 <- !I1
  
  if (any(I1)) {
    y[I1] <- (exp(z[I1]) - 1) / z[I1]
  }
  
  if (any(I2)) {
    zz <- z[I2]^2
    pv <- p[1] + zz * p[2]
    y[I2] <- 2 * pv / (q[1] + zz * (q[2] + zz * q[3]) - z[I2] * pv)
  }
  
  return(y)
}

# ---- tan2 function ----
tan2 <- function(x) {
  p <- c(129.221035, -8.87662377, 0.0528644456)
  q <- c(164.529332, -45.1320561, 1)
  pi4 <- pi / 4
  
  y <- numeric(length(x))
  xx <- abs(x)
  
  I1 <- xx > pi4
  I2 <- !I1
  
  if (any(I1)) {
    y[I1] <- tan(x[I1]) / x[I1]
  }
  
  if (any(I2)) {
    x2 <- (xx[I2] / pi4)^2
    y[I2] <- (p[1] + x2 * (p[2] + x2 * p[3])) / 
      (pi4 * (q[1] + x2 * (q[2] + x2 * q[3])))
  }
  
  return(y)
}

# ---- rstab function ----
rstab <- function(alpha, brime, n) {
  pi2 <- pi / 2
  thr1 <- 0.99
  eps <- 1 - alpha
  
  if (eps > -0.99) {
    tau <- brime / (pi2 * tan2(eps * pi2))
  } else {
    tau <- pi2 * alpha
    tau <- brime * eps * tau * tan2(tau)
  }
  
  a <- runif(n)
  phiby2 <- (a - 0.5) * pi2
  a <- phiby2 * tan2(phiby2)
  b <- phiby2 * eps
  bb <- tan2(b)
  b <- b * tan2(b)
  
  da <- a^2
  a2 <- 1 - da
  a2p <- 1 + da
  
  db <- b^2
  b2 <- 1 - db
  b2p <- 1 + db
  
  z <- -log(runif(n))
  z <- z * a2 * b2p
  z <- log(a2p * (b2 + phiby2 * bb * (2 * tau)) / z)
  z <- z / (1 - eps)
  d <- z * d2(z * eps)
  
  y <- a * b
  y <- 1 + y
  stab <- (1 + d * eps)
  stab <- (stab * 2) *
    ((a - b) * y - (phiby2 * tau) * bb * (b * a2 - a * 2)) /
    (a2 * b2p) + d * tau
  
  return(stab)
}

rsub <- function(alpha,n){
  # CHECK!
  # Simulates a standard stable subordinator variable X with
  # Laplace transform E[e^{-tX}] = exp{- t^alpha} which is also such that
  #
  # t*P(X>t^1/alpha) ~ 1/Gamma(1-alpha), as t -->\infty.
  #
  x = rstab(alpha,1,n);
  return((x-min(x))*(cos(alpha*pi/2)^(1/alpha)))
}

# Example usage:
# set.seed(42)
# samples <- rstab(alpha = 1.5, brime = 0, n = 1000)
# hist(samples, breaks = 50, main = "Alpha-stable random variates", col = "skyblue")
