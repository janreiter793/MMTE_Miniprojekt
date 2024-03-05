library(Matrix)
library(MatrixExtra)
library(magrittr)
set.seed(1)

# Profile likelihood procedure
profile.likelihood <- function(a, y, X, maximize = T) {
  if(maximize) {
    print(paste("Running profile likelihood with a =", round(a, 3)))
  }
    
  n <- length(y)
  
  # Obtain B^-1 on the basis of the given a
  i <- numeric(n * (n + 1) / 2)
  j <- numeric(n * (n + 1) / 2)
  x <- numeric(n * (n + 1) / 2)
  abs <- a^(1:n - 1)
  for(k in 1:n) {
    i[(((k - 1) * k / 2) + 1):(k * (k + 1) / 2)] <- rep(k, k)
    j[(((k - 1) * k / 2) + 1):(k * (k + 1) / 2)] <- 1:k
    x[(((k - 1) * k / 2) + 1):(k * (k + 1) / 2)] <- abs[1:k] %>% rev
  }
  B <- sparseMatrix(i = i, j = j, x = x, dims = c(n, n))
  Binv <- B %>% solve

  # Obtain sqrt(D^-1). D is diagonal, hence, symmetric. Therefore sqrt-matrix is
  # just the sqrt of all the entries of D.
  D <- sparseMatrix(i = 1:n, j = 1:n, 
                    x = c(1 / (1 - a^2), rep(1, n - 1)), 
                    dims = c(n, n))
  sqrtDinv <- D %>% sqrt %>% solve
  
  #Compute S
  S <- sqrtDinv %*% Binv
  
  ytilde <- as.numeric(S %*% y) # Some conversions of formats needed so that lm() is happy 
                                # (wants data to be of type numeric and design matrix Xtilde 
                                # to be of ordinary matrix type)
  
  Xtilde <- as.matrix(S %*% cbind(rep(1, n), X)) # Why add column of ones ? To include intercept 3 (beta_0)

  # Obtain beta and tau estimates
  fit <- lm(ytilde ~ -1 + Xtilde)
  
  #compute determinant of S to find log-likelihood in terms of Y.
  detS <- S %>% det
  logLikY <- logLik(fit) + log(detS)
  
  if (maximize) {
    #return likelihood of data y given a (NB log likelihood for ytilde can be extracted using logLik(fit))
    logLikY %>% as.numeric %>% return
  } else {
      #return likelihood of data y given a as well as fitted coeffiecients and variance
    return(list(log.likelihood = as.numeric(logLikY), 
                a              = a,
                beta           = coef(fit), 
                tau.sq         = sigma(fit)^2))
  }
}

#simulate data
n        <- 1000
a        <- 0.5
tau_sqrd <- 1

# Obtain B^-1 on the basis of the true parameter a
i <- numeric(n * (n + 1) / 2)
j <- numeric(n * (n + 1) / 2)
x <- numeric(n * (n + 1) / 2)
abs <- a^(1:n - 1)
for(k in 1:n) {
  i[(((k - 1) * k / 2) + 1):(k * (k + 1) / 2)] <- rep(k, k)
  j[(((k - 1) * k / 2) + 1):(k * (k + 1) / 2)] <- 1:k
  x[(((k - 1) * k / 2) + 1):(k * (k + 1) / 2)] <- abs[1:k] %>% rev
}
B <- sparseMatrix(i = i, j = j, x = x, dims = c(n, n))
Binv <- B %>% solve

# Obtain D^-1/2
D <- sparseMatrix(i = 1:n, j = 1:n, 
                  x = c(1 / (1 - a^2), rep(1, n - 1)), 
                  dims = c(n, n))
sqrtD <- D %>% sqrt

nu <- sqrt(tau_sqrd) * sqrtD %*% rnorm(n) # Droot inverse of Dinvroot and Binv defined as above
U <- solve(Binv, nu)                      # This corresponds to computing U=B eps

plot(U, type="l") ;grid() # Take a look at simulated errors
mean(U)
var(as.numeric(U))
acf(as.numeric(U)); grid()

# Generate data
x <- rnorm(n)
X <- matrix(x, ncol = 1)
y <- 3 + 2 * x + U

afit <- optimize(profile.likelihood, 
                 interval = c(-1, 1), 
                 y = y, X = X, 
                 maximum = T)

# Log-likelihood and ML estimates for, a, beta, and tau^2
profile.likelihood(afit$maximum, y, X, maximize = F)
