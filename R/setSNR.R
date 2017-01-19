setSNR <- function(noise, A){
  p <- ncol(A)
  inv <- solve(diag(p) - A)
  X <- noise%*%inv
  cvec <- numeric(p)
  for (j in 2:p) {
    ij <- 1:(j - 1)
    
    getLoss <- function(c){
      signalC <- X[, ij, drop = FALSE] %*% (c*A[ij, j, drop = FALSE])
      (1 - var(signalC + noise[,j]))^2
    }
    
    continue <- TRUE
    range <- sort(c(exp(seq(-2, log(50), length.out = 60)), seq(0.75, 1.25, by = 0.01)))
    range <- range[range != 0] 
    lossOld <- getLoss(1)
    
    while(continue){
      lossVals <- sapply(range, getLoss)
      c <- range[which.min(lossVals)]
      lossNew <- getLoss(c)
      # print(lossNew)
      # print(A)
      if(lossNew < lossOld){
        cvec[j] <- c
        A[ij, j] <- c%*%A[ij, j]
      }
      # print(A)
      # loss with new A
      loss <- getLoss(1)
      if(loss < 0.001 | (lossOld - lossNew) < 0.0001){
        continue <- FALSE
      }
      lossOld <- loss
      # print(loss)
    }
    
    X[, j] <- X[, ij, drop = FALSE] %*% A[ij, j] + noise[,j]
    
  }
  A
}