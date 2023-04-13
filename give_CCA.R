#----------------------
# Inputs:
#----------------------
## l = a list with the first item the initial alpha estimate and the second item the initial beta estimate
## C = the C parameter for the Hessian matrix
## X = the first set of variables
## Y = the second set of variables
## elements = which of the (p + q) parameters should the variance be calculated for
## nlC = the nodewise Lasso tuning parameter
#----------------------
# Outputs a list containing:
#----------------------
## matx = a data frame with debiased and initial x vectors in first and second columns, respectively
## maty = analagous data frame as matx but for the y vectors
## cond = the condition number of the Hessian matrix
## var_calc = a (p + q)-length vector with the variance of the parameters specified by the "elements" argument. Non-specified positions return "NA"
## rhosq_db = the debiased estimate of the squared first canonical correlation
## rhosq_var = the variance estimate of the debiased squared first canonical correlation

give_CCA <- function(alpha, beta, C, X, Y, elements, nlC)
{
 # a <- l[[1]] %>% as("sparseMatrix")
  #b <- l[[2]] %>% as("sparseMatrix")
  
  a <- alpha
  b <- beta
 
  Sx <- var(X)
  Sy <- var(Y)
  Sxy <- cov(X,Y)
  Syx <- t(Sxy)
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(X)
  if(missing(elements)) elements <- 1:(p+q)
  if(missing(nlC)) nlC <- sqrt(log(p+q)/n)
  
  # Getting an initial rho estimate 
  rho_init <- abs(sum(a*(Sxy%*%b)))
  
  xhat <- sqrt(rho_init)*a
  yhat <- sqrt(rho_init)*b
  xmat <- X
  ymat <- Y
  
  #Hessian
  H <- hessian(Sxy, Sx, Sy, a, b, rho_init, C = C) 
  cond <- 1/suppressWarnings(rcond(H)) 
    
  #Setting lambda
  lam = nlC * sqrt(log(p+q)/n) 
    
  #Nodewise Lasso
  list.vec <- parallel::mclapply(1:ncol(H), nodewise, lam=lam, H=H)
    
  #The estimate of the precision matrix (inverse of Hessian)
  Phi <- matrix(unlist(list.vec), byrow=TRUE, ncol= p+q)
    
  # Calculate the pseudovariables Z_j: j = 1...n (for variance calculation)
  Z <- rep(NA, n)
  var_calc <- rep(NA, p + q)
  
  for(i in elements){
    for (j in 1:n) {
      piece1 <- ((rho_init*t(Phi[(1:p), i]) + ((t(Phi[(1:p), i])%*% Sx %*% xhat)%*%t(xhat))) %*% (outer(as.numeric(xmat[j,]),as.numeric(xmat[j,]),"*")%*%xhat)) 
      piece2 <- t(Phi[(1:p), i])%*% as.numeric(xmat[j,]) %*% t(as.numeric(ymat[j,]))%*%yhat
      piece3 <- ((rho_init*t(Phi[(p+1):(p+q), i]) + ((t(Phi[(p+1):(p+q), i])%*% Sy %*% yhat)%*%t(yhat))) %*% (outer(as.numeric(ymat[j,]), as.numeric(ymat[j,]),"*") %*% yhat)) 
      piece4 <- t(xhat)%*% as.numeric(xmat[j,]) %*% t(as.numeric(ymat[j,])) %*% (Phi[(p+1):(p+q), i])
      Z[j] <- piece1 - piece2 + piece3 - piece4
    }
    var_calc[i] <- var(Z)
  } 
  
  #Gradient
  grad1 <- (2*rho_init*(Sx %*% xhat)) - 2*(Sxy %*% yhat)
  grad2 <- (2*rho_init*(Sy %*% yhat)) - 2*(Syx %*% xhat)
  temp <- c(grad1,grad2)
  print(dim(temp))
  # De-biasing
  #print(dim(cbind(xhat, yhat)))
  print(dim(Phi))
  de.bias <- c(xhat,yhat) - Phi%*%temp
  xhatdebias <- de.bias[1:p]
  yhatdebias <- de.bias[(p+1):(p+q)]
    
  matx <- data.frame(d.a=xhatdebias, i.a=as.vector(xhat))
  maty <- data.frame(d.b=yhatdebias, i.b=as.vector(yhat))
    
  xhat <- as.vector(xhat)
  yhat <- as.vector(yhat)
    
  #Calculating the de-biased squared first canonical correlation
  rhosq_db <- (t(xhat) %*% Sxy %*% yhatdebias) + (t(xhatdebias) %*% Sxy %*% yhat) - (t(xhat) %*% Sxy %*% yhat)
    
  #Calculating the variance of the above estimate
  rhosq_var <- vector(length = n)
  for (i in 1:n) {
      X_i <- as.numeric(xmat[i,])
      Y_i <- as.numeric(ymat[i,])
      rhosq_var[i] <- ((crossprod(xhat, X_i)^2) * (crossprod(yhat, Y_i)^2)) / n
  }
  
  rhosq_var <- sum(rhosq_var) - ((rho_init)^4)
  
  #Putting the output list together
  listestim <- list(matx, maty, cond, var_calc, rhosq_db, rhosq_var)
  return(listestim)
}
