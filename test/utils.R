# Compute the variance explained
variance_explained <- function(C, U){
  fve = 0
  k = dim(U)[2]
  for (i in 1:k){
    fve = fve + U[,i] %*% C %*% U[,i]
  }
  fve
}
# Compute the fraction of variance explained (variance explained normalized by the trace of the covariance/correlation matrix)
fraction_variance_explained <- function(C, U){
  fve = variance_explained(C,U)
  fve = fve / sum(diag(C))
  fve
}
# Compute the orthogonality violation
orthogonality_violation <- function(U){
  sum(abs(t(U) %*% U - diag(dim(U)[2])))
}
# Compress dataset
compress_data <- function(df, U){
  data.frame(as.matrix(df) %*% U)
}
