library(Matrix)
library(magrittr)

#reparametrization so we are sure to have positive phi and |a|<1
#theta=(logit((a + 1) / 2), log(phi))

#if A = t(R) %*% R then we can compute z=A^{-1}x by solving Az=x wrt. z. That is
#first temp=solve(t(R),x) and secondly z=solve(R,temp)

neg.log.profile.likelihood <- function(theta, y, X, minimize=T){
  a   <- 2 * exp(theta[1]) / (1 + exp(theta[1])) - 1
  phi <- exp(theta[2])
  if(minimize) print(c(a, phi))
  
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
  Winvy <- solve(CholQtilde, u)
  
  U <- solve(L, Q %*% X)
  WinvX <- solve(CholQtilde, U)
  
  betahat <- solve(t(X) %*% WinvX) %*% t(X) %*% Winvy
  
  #check when n small
  #Winv <- solve(Qtilde) %*% Q
  #betahat=solve(t(X)%*%Winv%*%X)%*%t(X)%*%Winv%*%y#OK
  
  #compute \hat \sigma^2(\theta)
  residual <- y - X %*% betahat
  u <- solve(L, Q %*% residual)
  WinvRes <- solve(CholQtilde, u)
  sigma2hat <- sum(residual * WinvRes) / n
  
  #check
  #z3 <- Winv %*% residual
  #sigma2hat <- sum(residual * z3) / n
  #sigma2hat
  
  #compute log likelihood. Note: determinant returns log determinant of original matrix.
  #check (only for small n)
  #V=sigma2hat*(diag(rep(1,n))+solve(Q)*phi)
  #-log(det(V))/2#-1.92
  #Vinv=solve(Qtilde)%*%Q/sigma2hat#OK
  
  #detVinv=det(Q)/(det(Qtilde)*sigma2hat^n)
  #logdetVinvhalf=log(det(Q))/2-log(det(Qtilde))/2-n*log(sigma2hat)/2
  
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

# Simulate data
n <- 1000000
a <- 0.5
tau2 <- 0.25
sigma2 <- 1
ii <- c(1:n, 2:n)
jj <- c(1:n, 1:(n - 1))
Binvij <- c(rep(1, n), rep(-a, n - 1))
Binv <- sparseMatrix(i = ii, j = jj, x = Binvij, dims = c(n, n))
Droot <- sparseMatrix(i = c(1:n), j = c(1:n), x = c(sqrt(1 / (1 - a^2)), rep(1, n - 1)))
nu <- sqrt(tau2) * Droot %*% rnorm(n)
U <- solve(Binv, nu)

x <- rnorm(n)
X <- cbind(rep(1, n), x)
y <- 3 + 2 * x + U + rnorm(n, sd = sqrt(sigma2))

tempa <- (a + 1) / 2
thetastart <- c(log(tempa / (1 - tempa)), log(tau2 / sigma2))
neg.log.profile.likelihood(thetastart, y, X, minimize = F)

tempa <- (0.3 + 1) / 2
thetastart <- c(log(tempa / (1 - tempa)), log(1))
thetafit <- optim(thetastart, neg.log.profile.likelihood, y = y, X = X)
thetafit
neg.log.profile.likelihood(c(thetafit$par[1], thetafit$par[2]), y, X, minimize = F)

2 * exp(thetafit$par[1]) / (1 + exp(thetafit$par[1])) - 1
exp(thetafit$par[2])

