library(mvtnorm)
#Simulate standard normal data matrix: first generate alpha and beta
p <- 50; q <- 50; al <- c(rep(1, 10), rep(0, 40));
be <- c(rep(0,25), rnorm(25,1))
#Normalize alpha and beta
al <- al/sqrt(sum(al^2))
be <- be/sqrt(sum(be^2))
n <- 300; rho <- 0.5
#Creating the covariance matrix
Sigma_mat <- function(p,q,al,be, rho)
{
  Sx <- diag(rep(1,p), p, p)
  Sy <- diag(rep(1,q), q, q)
  Sxy <- tcrossprod(crossprod(rho*Sx, outer(al, be)), Sy)
  Syx <- t(Sxy)
  rbind(cbind(Sx, Sxy), cbind(Syx, Sy))
}
truesigma <- Sigma_mat(p,q,al,be, rho)
#Simulating the data
Z <- mvtnorm::rmvnorm(n, sigma = truesigma)
x <- Z[,1:p]
y <- Z[,(p+1):(p+q)]
elements <- 1:p
nlC <- log(p+q)/n