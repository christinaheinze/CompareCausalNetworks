runBivariateANM<- function(X, parentsOf, variableSelMat, pointEst, verbose, 
                           result, ...){
  
  dots <- list(...)
  if(length(dots) > 0){
    warning("options provided via '...' not taken")
  }
  
  bivanmmat <- bivariateANM(X, parentsOf =parentsOf, 
                            variableSelMat = variableSelMat, 
                            silent = !verbose)
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- (wh <- which(bivanmmat$causalParents[, parentsOf[k]]>0))
    if(pointEst) 
      attr(result[[k]],"coefficients") <- bivanmmat$scoreMat[ wh,parentsOf[k] ]
  }
  
  result
}