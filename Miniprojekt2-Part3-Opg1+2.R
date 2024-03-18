library(Matrix)
library(MatrixExtra)
library(magrittr)

# Parameters
N_SIM <- 100

# Profile likelihood procedure for sigma = 0
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
  
  Xtilde <- as.matrix(S %*% cbind(rep(1, n), as.matrix(X))) # Why add column of ones ? To include intercept 3 (beta_0)
  
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

# Profile log likelihood for sigma != 0
neg.log.profile.likelihood <- function(theta, y, X, minimize=T){
  a   <- 2 * exp(theta[1]) / (1 + exp(theta[1])) - 1
  phi <- exp(theta[2])
  
  n <- length(y)
  
  # Construct Binv and Dinv as sparse matrices
  ii <- c(1:n, 2:n)
  jj <- c(1:n, 1:(n - 1))
  Binvij <- c(rep(1, n), rep(-a, n - 1))
  Binv <- sparseMatrix(i = ii, j = jj, x = Binvij, dims = c(n, n))
  Dinv <- sparseMatrix(i = c(1:n), j = c(1:n), x = c(1 - a^2, rep(1, n - 1)))
  
  #Compute Q and Qtilde (again as sparse matrices
  Q <- t(Binv) %*% Dinv %*% Binv
  Qtilde <- sparseMatrix(i = 1:n, j = 1:n, x = rep(phi, n), dims = c(n, n)) + Q
  
  # compute Cholesky factorization
  CholQtilde <- chol(Qtilde)
  L <- t(CholQtilde)
  
  #compute \hat \beta(\theta)
  #(X^T W^{-1} X)^{-1} X^\T W^{-1} y
  #using Cholesky factor CholQtilde
  u <- solve(L, Q %*% y)
  Winvy <- solve(CholQtilde) %*% u
  
  Xtilde <- as.matrix(cbind(rep(1, n), as.matrix(X)))
  
  U <- solve(L, Q %*% Xtilde)
  WinvX <- solve(CholQtilde) %*% U
  
  betahat <- solve(t(Xtilde) %*% WinvX) %*% t(Xtilde) %*% Winvy
  
  #compute \hat \sigma^2(\theta)
  residual <- y - Xtilde %*% betahat
  u <- solve(L, Q %*% residual)
  WinvRes <- solve(CholQtilde) %*% u
  sigma2hat <- sum(residual * WinvRes) / n
  
  detQ <- det(Q)
  loglikelihood <- -n * log(sigma2hat) / 2 - determinant(CholQtilde)$modulus + 
    log(detQ) / 2 - (2 * sigma2hat)^(-1) * sum(residual * WinvRes)
  
  loglikelihood %<>% as.numeric
  # return negative log likelihood for later use with optim() which minimizes as default.
  if (minimize)
    return(-loglikelihood)
  else
    return(list(loglikelihood, betahat, sigma2hat))
}

logLik_restrictedModel <- function(y, X, show.estimates = F) {
  afit <- optimize(profile.likelihood, interval = c(-1, 1), 
                   y = y, X = X, maximum = T)
  profile.likelihood(afit$maximum, y, X, maximize = (!show.estimates)) %>% return
}

logLik_unrestrictedModel <- function(y, X, show.estimates = F) {
  thetafit <- optim(c(0.5, 0), neg.log.profile.likelihood, y = y, X = X)
  if(show.estimates) {
    neg.log.profile.likelihood(c(thetafit$par[1], thetafit$par[2]), y, X, minimize = FALSE) %>% return
  } else {
    -neg.log.profile.likelihood(c(thetafit$par[1], thetafit$par[2]), y, X) %>% return
  }
}

# Load data table
data <- read.table("C:\\Users\\janre\\Documents\\uni\\8. Semester\\MMTE\\Miniprojekt_2\\consumption.txt", 
                   sep = "" , header = T)
data %<>% as.data.frame 
data$D %<>% factor
data

y <- data$elforbrug
X <- data[,c(2,3)] %>% as.data.frame

restricted <- logLik_restrictedModel(y, X, show.estimates = TRUE)
unrestricted <- logLik_unrestrictedModel(y, X, show.estimates = TRUE)

Q <- -2 * (restricted$log.likelihood - unrestricted[[1]])

beta_hat <- restricted$beta %>% unname
a_hat    <-restricted$a
tau2_hat <- restricted$tau.sq
n <- y %>% length

mu <- beta_hat[1] * rep(1, n) + beta_hat[2] * X$temperatur + beta_hat[3] * as.numeric(X$D)
ii <- c(1:n, 2:n)
jj <- c(1:n, 1:(n - 1))
Binvij <- c(rep(1, n), rep(-a_hat, n - 1))
Binv <- sparseMatrix(i = ii, j = jj, x = Binvij, dims = c(n, n))

# Obtain D^-1/2
D <- sparseMatrix(i = 1:n, j = 1:n, 
                  x = c(1 / (1 - a_hat^2), rep(1, n - 1)), 
                  dims = c(n, n))
sqrtD <- D %>% sqrt

Qs <- numeric(N_SIM)
for(i in 1:N_SIM) {
  print(paste("Running simulation no.", i, "out of", N_SIM))
  nu <- sqrt(tau2_hat) * sqrtD %*% rnorm(n) # Droot inverse of Dinvroot and Binv defined as above
  U <- solve(Binv, nu)                      # This corresponds to computing U=B eps
  y_simulated <- mu + U
  
  restricted_sim <- logLik_restrictedModel(y_simulated, X)
  unrestricted_sim <- logLik_unrestrictedModel(y_simulated, X)
  Qs[i] <- -2 * (restricted_sim - unrestricted_sim)
}

p_val <- 1 - mean(Qs < Q)
