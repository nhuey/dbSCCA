# In this script, I'm going to lay out the flow for the example file for the paper code

library(mvtnorm)
file.sources = list.files("~/Dropbox/Sparse CCA/Annotated Codes/Auxillary_Functions_from_Nathan/",pattern="*.R", full.names = TRUE)
sapply(file.sources[-2],source,.GlobalEnv)
#Simulate standard normal data matrix: first generate alpha and beta
p <- 500; q <- 100; s1 = 50; s2 = 10
set.seed(072024)
al <- rep(0, p); be <- rep(0,q)
#alpha.non.zero.index <- sample(1:p, size = s1, replace = FALSE)
alpha.non.zero.index <- 1:45

al[alpha.non.zero.index] <- rnorm(s1)
beta.non.zero.index <- sample(1:q, size = s2, replace = FALSE)
be[beta.non.zero.index] <- rnorm(s2)
#Normalize alpha and beta
al <- al/sqrt(sum(al^2))
be <- be/sqrt(sum(be^2))
n <- 10000; rho <- 0.7
#Creating the covariance matrix

# Maybe we got good AUC even with irrelevant CNAs because of the correlation within the CNAs. This can be tested. In our fake data, can you change the covariance matrix of X so that everyone has some correlation (say 0.2, 0.4, 0.6, etc.). Then rerun everything with the wrong version of mRMR and check the AUC. How much AUC do we get and how does it change with correlation?
Sigma_mat <- function(p,q,al,be, rho)
{
  Sx <- diag(rep(1-0.1 ,p), p, p) + matrix(0.1, nrow = p, ncol = p)
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
elements <- 1:(p+q)
nlC <- log(p+q)/n

# First step: MRMR
# Choose 45 non-zero CNA's (X) to be involved in ER-status and 5 zero CNA's
non.zero.CNA <- sample(alpha.non.zero.index, 45, replace = FALSE)
zero.CNA <- sample((1:p)[-alpha.non.zero.index], size = 1, replace = FALSE)

er.coeff <- rep(0, p)
er.coeff[c(non.zero.CNA, zero.CNA)] <- 1

z = 1 + x %*% er.coeff
pr = 1/(1+exp(-z))
er = rbinom(n,1,pr)

data.x <- as.data.frame(x)
data.x$er <- er
#logis.reg <- glm(er~., data = data.x)
#summary((logis.reg))
# ER-status is now simulated depending only on the CNA
#sort(c(non.zero.CNA, zero.CNA))
#sort(mrmr.selection.final)
# Choose CNAs based on MRMR analysis of ER-status
CNA.indices <- choose.CNAs(x, er, size = 70)
#sort(CNA.indices)

# Compare AUCs of current method and old method - maybe do 3 splits?

mai.scca <- SCCA(x[, CNA.indices], y, lambda.alpha = 0.1, lambda.beta = 0.1)
example.result <- give_CCA(alpha = mai.scca$beta, beta = mai.scca$alpha, C = 2, X = x[,CNA.indices], Y = y, elements = 1:170, nlC)
example.result[[6]]



x6 <- (example.result)[[1]][,1]
x6.var <- (example.result)[[4]][1:70]
z6 <- x6/sqrt(4*x6.var/n)
p.val.6 <- 2*pnorm(abs(z6), lower.tail = FALSE)

index.6 <-which(p.adjust(p.val.6) < 0.05)
sort(CNA.indices[index.6])


# our method didn't select any false signals and only missed 10 out of 45
sort(CNA.indices)
# actual non-zero: 1-45 (X) and 
sort(beta.non.zero.index)
#index.6 <- which(p.val.6 < 2.5e-6)

# Just going w/ SCCA
sparse.x <- (example.result)[[1]][,2]
sparse.index <- which(sparse.x != 0)
sort(CNA.indices[sparse.index]) # just going with this misses 33 more CNA

y6 <- (example.result)[[2]][,1]
y6.var <- (example.result)[[4]][71:170]
zy6 <- y6/sqrt(4*y6.var/n)
p.val.6.y <- 2*pnorm(abs(zy6), lower.tail = FALSE)
#index.6.y <- which(p.val.6.y < 2.5e-6)

index.6.y <-which(p.adjust(p.val.6.y) < 0.05)
sort(index.6.y)
###
sparse.y <- (example.result)[[2]][,2]
sparse.index.y <- which(sparse.y != 0)
sort(sparse.index.y) # just going with this misses one of the genes


# Wrap into single function


cna.data

old.indices <- read.table("~/Dropbox/Breast_cancer/CNA_selection/CNA_final_selections/mrmr_final_select_fulldata.txt")

new.indices <- read.table("~/Dropbox/Breast_cancer/CNA_selection/CNA_final_selections/mrmr_final_select_fulldata_new.txt")


cna.corrs <- cor(cna.data[,old.indices$x], cna.data[,new.indices$x])

max.corrs <- apply(cna.corrs, 2, max, na.rm = TRUE )


histogram(c(cna.corrs))

min(max.corrs[!is.infinite(max.corrs)])
sum(abs(max.corrs) > 0.8)

cna.corrs[2,1]


