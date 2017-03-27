runBivariateCAM<- function(X, parentsOf, variableSelMat, pointEst, verbose, 
                           ...){
  
  dots <- list(...)
  if(length(dots) > 0){
    warning("options provided via '...' not taken")
  }
  
  bivcammat <- bivariateCAM(X, parentsOf=parentsOf, 
                            variableSelMat=variableSelMat, 
                            silent = !verbose)
  
  result <- vector("list", length = length(parentsOf))
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- (wh <- which(bivcammat$causalParents[, parentsOf[k]]>0))
    
    attr(result[[k]],"parentsOf") <- parentsOf[k]
    
    if(pointEst)
      attr(result[[k]],"coefficients") <- bivcammat$scoreMat[ wh,parentsOf[k] ]
  }
  
  resMat <- bivcammat$causalParents
  if(pointEst)
    resMat[resMat == TRUE] <- bivcammat$scoreMat[resMat == TRUE]
  
  list(resList = result, resMat = resMat)
}