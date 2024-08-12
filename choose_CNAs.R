library(mRMRe)
choose.CNAs <- function(data, outcomes, size = ncol(data)){
cna.data <- data
# Number of desired features
k <- size

# mRMR
df <- data.frame(cbind(outcomes, cna.data))
# converting columns from integer to numeric type
for(i in 1:ncol(df)){df[,i] <- as.numeric(df[,i])}

dd <- mRMR.data(data = df)
mrmr <- mRMR.ensemble(data = dd, target_indices = 1, solution_count = 40,
                      feature_count = ncol(cna.data)-1)
# Subtract b/c we added the outcome into the data frame
# mrmr.selection <- (((mrmr@filters)[[1]]))-1
# mrmr.selection.ranks <- apply(mrmr@filters[[1]]-1, 2, function(x) match(1:(length(x)+1), rev(x)))
# mrmr.selection.ranks[is.na(mrmr.selection.ranks)] <- nrow(mrmr.selection.ranks) + 1
# mrmr.selection.rank <- apply(mrmr.selection.ranks, 1, mean)
# mrmr.selection.final <- match(sort(mrmr.selection.rank)[1:k], mrmr.selection.rank)
#browser()
mrmr.selection <- solutions(mrmr)$'1'
mrmr.selection.ranks <- apply(mrmr.selection, 2, function(x) match(1:(length(x)+1), x))
#mrmr.selection.ranks <- apply(mrmr@filters[[1]]-1, 2, function(x) match(1:(length(x)+1), x))
mrmr.selection.ranks[is.na(mrmr.selection.ranks)] <- nrow(mrmr.selection.ranks)
mrmr.selection.rank <- apply(mrmr.selection.ranks, 1, mean)
mrmr.selection.final <- match(sort(mrmr.selection.rank)[1:k], mrmr.selection.rank)
# I believe we need to adjust the column numbers here
mrmr.selection.final <- mrmr.selection.final - 1

return(mrmr.selection.final)
}

choose.CNAs.old <- function(data, outcomes, size = ncol(data)){
  cna.data <- data
  # Number of desired features
  k <- size
  
  # mRMR
  df <- data.frame(cbind(outcomes, cna.data))
  # converting columns from integer to numeric type
  for(i in 1:ncol(df)){df[,i] <- as.numeric(df[,i])}
  
  dd <- mRMR.data(data = df)
  mrmr <- mRMR.ensemble(data = dd, target_indices = c(1), solution_count = 40,
                        feature_count = ncol(cna.data)-1)
  # Subtract b/c we added the outcome into the data frame
  mrmr.selection <- (((mrmr@filters)[[1]]))
  mrmr.selection.ranks <- apply(mrmr@filters[[1]]-1, 2, function(x) match(1:(length(x)+1), x))
  mrmr.selection.ranks[is.na(mrmr.selection.ranks)] <- nrow(mrmr.selection.ranks) + 1
  mrmr.selection.rank <- apply(mrmr.selection.ranks, 1, mean)
  mrmr.selection.final <- match(sort(mrmr.selection.rank)[1:k], mrmr.selection.rank)
  
  
  mrmr.selection[,1]-1
  
  return(mrmr.selection.final)
}
