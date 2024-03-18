library(Matrix)
library(MatrixExtra)
library(magrittr)
library(future.apply)
library(microbenchmark)
set.seed(1)

# Parameters
TAU2  <- 1
BETA  <- c(3, 2)
SIM_N <- 1000     # Takes roughly 1½ minutes to run with SIM_N <- 1000
BENCHMARK <- FALSE
BENCHMARK_ITER <- 1

# Profile likelihood procedure
profile.likelihood <- function(a_hat, y, X, maximize = T) {
  n <- length(y)
  
  # Obtain B^-1 on the basis of the given a
  ii <- c(1:n, 2:n)
  jj <- c(1:n, 1:(n - 1))
  Binvij <- c(rep(1, n), rep(-a_hat, n - 1))
  Binv <- sparseMatrix(i = ii, j = jj, x = Binvij, dims = c(n, n))
  
  # Obtain sqrt(D^-1). D is diagonal, hence, symmetric. Therefore sqrt-matrix is
  # just the sqrt of all the entries of D.
  D <- sparseMatrix(i = 1:n, j = 1:n, 
                    x = c(1 / (1 - a_hat^2), rep(1, n - 1)), 
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
                a              = a_hat,
                beta           = coef(fit), 
                tau.sq         = sigma(fit)^2))
  }
}

MLestimate <- function(y, X) {
  afit <- optimize(profile.likelihood, interval = c(-1, 1), 
                   y = y, X = X, maximum = T)
  profile.likelihood(afit$maximum, y, X, maximize = F) %>% return
}

simulate <- function(a, n, beta = BETA, tau2 = TAU2, N = SIM_N) {
  a_estimates     <- numeric(N)
  beta1_estimates <- numeric(N)
  beta2_estimates <- numeric(N)
  tau2_estimates  <- numeric(N)
  
  # Obtain B^-1 on the basis of the true parameter a
  ii <- c(1:n, 2:n)
  jj <- c(1:n, 1:(n - 1))
  Binvij <- c(rep(1, n), rep(-a, n - 1))
  Binv <- sparseMatrix(i = ii, j = jj, x = Binvij, dims = c(n, n))
  
  # Obtain D^-1/2
  D <- sparseMatrix(i = 1:n, j = 1:n, 
                    x = c(1 / (1 - a^2), rep(1, n - 1)), 
                    dims = c(n, n))
  sqrtD <- D %>% sqrt
  
  for(i in 1:N) {
    nu <- sqrt(tau2) * sqrtD %*% rnorm(n) # Droot inverse of Dinvroot and Binv defined as above
    U <- solve(Binv, nu)                  # This corresponds to computing U=B eps
    
    # Generate data
    x <- rnorm(n)
    X <- matrix(x, ncol = 1)
    y <- beta[1] + beta[2] * x + U
    
    res <- MLestimate(y, X)
    a_estimates[i]     <- res$a
    beta1_estimates[i] <- res$beta[1]
    beta2_estimates[i] <- res$beta[2]
    tau2_estimates[i]  <- res$tau.sq
  }
  
  return(list(a = a_estimates, beta1 = beta1_estimates, 
              beta2 = beta2_estimates, tau2 = tau2_estimates))
}

run_sim <- function(params) {
  return(simulate(params[1], params[2]))
} 

parameters_comb <- list(c(0, 20), c(0.5, 20), c(0.99, 20),
                        c(0, 1000), c(0.5, 1000), c(0.99, 1000))
plan(multisession, workers = 6)

if(!BENCHMARK) {
  res <- future_lapply(parameters_comb, run_sim) %>% suppressWarnings
} else if(BENCHMARK) {
  microbenchmark::microbenchmark(
    res <- future_lapply(parameters_comb, run_sim) %>% suppressWarnings,
    times = BENCHMARK_ITER
  )
}

par(mfrow = c(2, 2))
# Simulation study with a = 0, n = 20
res[[1]]$a %T>% qqnorm(main = "Distribution for a_hat, (a = 0, n = 20)") %>% qqline; grid()
res[[1]]$beta1 %T>% qqnorm(main = "Distribution for beta1_hat, (a = 0, n = 20)") %>% qqline; grid()
res[[1]]$beta2 %T>% qqnorm(main = "Distribution for beta2_hat, (a = 0, n = 20)") %>% qqline; grid()
res[[1]]$tau2 %T>% qqnorm(main = "Distribution for tau2_hat, (a = 0, n = 20)") %>% qqline; grid()

res[[1]]$a %>% hist(main = "Distribution for a_hat, (a = 0, n = 20)", breaks = 50); grid()
res[[1]]$beta1 %>% hist(main = "Distribution for beta1_hat, (a = 0, n = 20)", breaks = 50); grid()
res[[1]]$beta2 %>% hist(main = "Distribution for beta2_hat, (a = 0, n = 20)", breaks = 50); grid()
res[[1]]$tau2 %>% hist(main = "Distribution for tau2_hat, (a = 0, n = 20)", breaks = 50); grid()

# Simulation study with a = 0.5, n = 20
res[[2]]$a %T>% qqnorm(main = "Distribution for a_hat, (a = 0.5, n = 20)") %>% qqline; grid()
res[[2]]$beta1 %T>% qqnorm(main = "Distribution for beta1_hat, (a = 0.5, n = 20)") %>% qqline; grid()
res[[2]]$beta2 %T>% qqnorm(main = "Distribution for beta2_hat, (a = 0.5, n = 20)") %>% qqline; grid()
res[[2]]$tau2 %T>% qqnorm(main = "Distribution for tau2_hat, (a = 0.5, n = 20)") %>% qqline; grid()

res[[2]]$a %>% hist(main = "Distribution for a_hat, (a = 0.5, n = 20)", breaks = 50); grid()
res[[2]]$beta1 %>% hist(main = "Distribution for beta1_hat, (a = 0.5, n = 20)", breaks = 50); grid()
res[[2]]$beta2 %>% hist(main = "Distribution for beta2_hat, (a = 0.5, n = 20)", breaks = 50); grid()
res[[2]]$tau2 %>% hist(main = "Distribution for tau2_hat, (a = 0.5, n = 20)", breaks = 50); grid()

res[[2]]$a %>% mean

# Simulation study with a = 0.99, n = 20
res[[3]]$a %T>% qqnorm(main = "Distribution for a_hat, (a = 0.99, n = 20)") %>% qqline; grid()
res[[3]]$beta1 %T>% qqnorm(main = "Distribution for beta1_hat, (a = 0.99, n = 20)") %>% qqline; grid()
res[[3]]$beta2 %T>% qqnorm(main = "Distribution for beta2_hat, (a = 0.99, n = 20)") %>% qqline; grid()
res[[3]]$tau2 %T>% qqnorm(main = "Distribution for tau2_hat, (a = 0.99, n = 20)") %>% qqline; grid()

res[[3]]$a %>% hist(main = "Distribution for a_hat, (a = 0.99, n = 20)", breaks = 50); grid()
res[[3]]$beta1 %>% hist(main = "Distribution for beta1_hat, (a = 0.99, n = 20)", breaks = 50); grid()
res[[3]]$beta2 %>% hist(main = "Distribution for beta2_hat, (a = 0.99, n = 20)", breaks = 50); grid()
res[[3]]$tau2 %>% hist(main = "Distribution for tau2_hat, (a = 0.99, n = 20)", breaks = 50); grid()

# Simulation study with a = 0, n = 1000
res[[4]]$a %T>% qqnorm(main = "Distribution for a_hat, (a = 0, n = 1000)") %>% qqline; grid()
res[[4]]$beta1 %T>% qqnorm(main = "Distribution for beta1_hat, (a = 0, n = 1000)") %>% qqline; grid()
res[[4]]$beta2 %T>% qqnorm(main = "Distribution for beta2_hat, (a = 0, n = 1000)") %>% qqline; grid()
res[[4]]$tau2 %T>% qqnorm(main = "Distribution for tau2_hat, (a = 0, n = 1000)") %>% qqline; grid()

res[[4]]$a %>% hist(main = "Distribution for a_hat, (a = 0, n = 1000)", breaks = 50); grid()
res[[4]]$beta1 %>% hist(main = "Distribution for beta1_hat, (a = 0, n = 1000)", breaks = 50); grid()
res[[4]]$beta2 %>% hist(main = "Distribution for beta2_hat, (a = 0, n = 1000)", breaks = 50); grid()
res[[4]]$tau2 %>% hist(main = "Distribution for tau2_hat, (a = 0, n = 1000)", breaks = 50); grid()

# Simulation study with a = 0.5, n = 1000
res[[5]]$a %T>% qqnorm(main = "Distribution for a_hat, (a = 0.5, n = 1000)") %>% qqline; grid()
res[[5]]$beta1 %T>% qqnorm(main = "Distribution for beta1_hat, (a = 0.5, n = 1000)") %>% qqline; grid()
res[[5]]$beta2 %T>% qqnorm(main = "Distribution for beta2_hat, (a = 0.5, n = 1000)") %>% qqline; grid()
res[[5]]$tau2 %T>% qqnorm(main = "Distribution for tau2_hat, (a = 0.5, n = 1000)") %>% qqline; grid()

res[[5]]$a %>% hist(main = "Distribution for a_hat, (a = 0.5, n = 1000)", breaks = 50); grid()
res[[5]]$beta1 %>% hist(main = "Distribution for beta1_hat, (a = 0.5, n = 1000)", breaks = 50); grid()
res[[5]]$beta2 %>% hist(main = "Distribution for beta2_hat, (a = 0.5, n = 1000)", breaks = 50); grid()
res[[5]]$tau2 %>% hist(main = "Distribution for tau2_hat, (a = 0.5, n = 1000)", breaks = 50); grid()

# Simulation study with a = 0.99, n = 1000
res[[6]]$a %T>% qqnorm(main = "Distribution for a_hat, (a = 0.99, n = 1000)") %>% qqline; grid()
res[[6]]$beta1 %T>% qqnorm(main = "Distribution for beta1_hat, (a = 0.99, n = 1000)") %>% qqline; grid()
res[[6]]$beta2 %T>% qqnorm(main = "Distribution for beta2_hat, (a = 0.99, n = 1000)") %>% qqline; grid()
res[[6]]$tau2 %T>% qqnorm(main = "Distribution for tau2_hat, (a = 0.99, n = 1000)") %>% qqline; grid()

res[[6]]$a %>% hist(main = "Distribution for a_hat, (a = 0.99, n = 1000)", breaks = 50); grid()
res[[6]]$beta1 %>% hist(main = "Distribution for beta1_hat, (a = 0.99, n = 1000)", breaks = 50); grid()
res[[6]]$beta2 %>% hist(main = "Distribution for beta2_hat, (a = 0.99, n = 1000)", breaks = 50); grid()
res[[6]]$tau2 %>% hist(main = "Distribution for tau2_hat, (a = 0.99, n = 1000)", breaks = 50); grid()
