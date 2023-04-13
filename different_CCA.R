# Different types of CCA estimators are here

#--------------------------------- Daniela 's method ----------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

# BiocManager::install("impute")
library("PMA")
#Not good if Sigma is not identity
cca_dan <- function(xmat,ymat)
{
  dxmat <- scale(xmat)
  dymat <- scale(ymat)
  dcca <- PMA::CCA(dxmat, dymat, typex = "standard", typez = "standard", standardize = FALSE)
  da <- dcca$u
  db <- dcca$v
  drho <- sum(t(da)%*%cov(dxmat,dymat)%*%db)
  list(da,db,drho)
}

#dv <- dcca$v.init  # the first k factors of the SVD that Daniella asked to save if doing
#this multiple time, to be passed when we use this next time

#------- normal naive CCA--------------------------
cca_pati <- function(xmat, ymat)
{
  pcca <- cancor(xmat, ymat)
  pa <- pcca$xcoef[1,]
  pb <- pcca$ycoef[1,]
  prho <- pcca$cor[1]
  list(pa, pb, prho)
}


#------ SCCA (MAI 2017)---------------------------
cca.mai <- function(xmat, ymat)
{
  n <- nrow(xmat)
  p <- ncol(xmat)
  q <- ncol(ymat)
  s.cca <- SCCA(xmat, ymat, lambda.alpha=sqrt(log(p)/n), lambda.beta = sqrt(log(q)/n))
  sb <- s.cca$alpha
  sa <-s.cca$beta
  srho <- sum(t(sa)%*%cov(xmat,ymat)%*%sb)
  list(sa, sb, srho)
}

# This function allows the tuning parameters of Mai's method to be changed
cca.mai2 <- function(xmat, ymat, aC, bC)
{
  n <- nrow(xmat)
  p <- ncol(xmat)
  q <- ncol(ymat)
  s.cca <- SCCA(xmat, ymat, lambda.alpha= aC*sqrt(log(p)/n), lambda.beta = bC*sqrt(log(q)/n))
  sb <- s.cca$alpha
  sa <-s.cca$beta
  srho <- sum(t(sa)%*%cov(xmat,ymat)%*%sb)
  list(sa, sb, srho)
}


#--------- WCCA (Willms 201?)-----------------------
source("~/Dropbox/Sparse_CCA/Annotated Codes/SAR.R")
# Takes too long to run 
cca_sar <- function(xmat, ymat)
{
  wcca <- SparseCCA(xmat, ymat, rank = 1)
  wa <- wcca$ALPHA
  wb <- wcca$BETA
  wrho <- wcca$cancors
  lama <- wcca$lambdaA
  lamb <- wcca$lambdaB
  list(wa, wb)
}


#-------- Parkhomenko CCA ----------
cca_park <- function(xmat, ymat)
{
  pkcca <- SCCA_Parkhomenko(xmat, ymat)
  pka <- pkcca$a
  pkb <- pkcca$b
  pkrho <- pkcca$cancor
  list(pka, pkb, pkrho)
}


#------- Waaijenborg CCA ---------
cca_waa <- function(xmat, ymat)
{
  wjcca <- Waaijenborg(xmat, ymat, rank=1)
  wja <- wjcca$uhat
  wjb <- wjcca$vhat
  wjrho <- wjcca$cancors
}



