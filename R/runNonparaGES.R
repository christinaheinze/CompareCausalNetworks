runNonparanormalGES <- function(X, parentsOf, variableSelMat, setOptions, directed, verbose, 
                                  result, ...){
  
  given.cov.mat <- cov(X)
  p <- ncol(X)
  n <- nrow(X)
  e <- eigen(given.cov.mat)
  sqrt.given.cov.mat <- e$vectors%*%sqrt(diag(e$values))
  
  ## generate from N(0,I)
  dat <- matrix(rnorm(p*n), ncol = p, nrow = n)
  
  ## transform data so that cov(dat) = given.cov.mat
  samp.cov.mat <- 2*sin(cor(dat,method="spearman")*pi/6)
  e <- eigen(samp.cov.mat, symmetric = T)
  e$values[which(e$values<0)] <- 0
  samp.cov.mat <- e$vectors%*%diag(e$values)%*%t(e$vectors)
  e <- eigen(samp.cov.mat)
  sqrt.samp.cov.mat <- e$vectors%*%sqrt(diag(e$values))
  X <- t(sqrt.given.cov.mat%*%solve(sqrt.samp.cov.mat,t(dat)))
  
  runGES(X, parentsOf, variableSelMat, setOptions, directed, verbose, 
         result, ...)
}