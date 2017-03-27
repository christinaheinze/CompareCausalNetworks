runBivariateANM<- function(X, parentsOf, variableSelMat, pointEst, verbose, 
                           ...){
  
  dots <- list(...)
  if(length(dots) > 0){
    warning("options provided via '...' not taken")
  }
  
  bivanmmat <- bivariateANM(X, parentsOf =parentsOf, 
                            variableSelMat = variableSelMat, 
                            silent = !verbose)
  
  result <- vector("list", length = length(parentsOf))
  for (k in 1:length(parentsOf)){
    result[[k]] <- (wh <- which(bivanmmat$causalParents[, parentsOf[k]]>0))
    
    attr(result[[k]],"parentsOf") <- parentsOf[k]
    
    if(pointEst) 
      attr(result[[k]],"coefficients") <- bivanmmat$scoreMat[ wh,parentsOf[k] ]
  }
  
  resMat <- bivanmmat$causalParents
  if(pointEst)
    resMat[resMat == TRUE] <- bivanmmat$scoreMat[resMat == TRUE]
 
  list(resList = result, resMat = resMat)
}