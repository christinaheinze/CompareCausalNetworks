runLINGAM <- function(X, parentsOf, pointEst, setOptions, verbose, ...){
  
  dots <- list(...)
  if(length(dots) > 0){
    warning("options provided via '...' not taken")
  }
  
  if(nrow(X)<=ncol(X)) 
    stop("LINGAM not suitable for high-dimensional data; 
         need nrow(X) > ncol(X)")
  
  res <- pcalg::lingam(X, verbose=verbose)
  B <- res$Bpruned
  B[B != 0] <- 1
  lingammat <- t(B)
  lingammatCoef <- t(res$Bpruned) 

  result <- vector("list", length = length(parentsOf))
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- (wh <- which(lingammat[, parentsOf[k]] == 1))
    attr(result[[k]],"parentsOf") <- parentsOf[k]
    
    if(pointEst)
      attr(result[[k]],"coefficients") <- lingammatCoef[ wh,parentsOf[k]]
  }
  
  if(length(parentsOf) < ncol(X)){
    lingammat <- lingammat[,parentsOf]
    lingammatCoef <- lingammatCoef[,parentsOf]
  }
  
  list(resList = result, resMat = if(pointEst) lingammatCoef else lingammat)
}