runRegression <- function(X, parentsOf, variableSelMat, pointEst, setOptions, 
                          directed, verbose, ...){
  
  dots <- list(...)
  if(length(dots) > 0){
    warning("options provided via '...' not taken")
  }
  
  result <- vector("list", length = length(parentsOf))
  
  # additional options for regression
  optionsList <- list("selfselect"=NULL)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
 
  for (k in 1:length(parentsOf)){
    if(round(k/100)==(k/100)) cat(" ",k)
    
    possibleVar <- (1:ncol(X))
    
    possibleVar <- (1:ncol(X))
    removeVar <- parentsOf[k]
    if(!is.null(variableSelMat) & is.null(optionsList$selfselect)){
      if(ncol(variableSelMat)==length(parentsOf)) 
        selc <- k 
      else 
        selc <- which( (1:ncol(X)) == parentsOf[k])
      removeVar <- 
        unique(c(removeVar,which( !variableSelMat[,selc] )))
    }else{
      if(!is.null(optionsList$selfselect)){
        gl <- glmnet::glmnet(X[, possibleVar[-removeVar],drop=FALSE], 
                             as.numeric(X[, parentsOf[k]]))
        nnz <- apply(coef(gl)!=0, 2,sum)
        beta <- 
          coef(gl, s=gl$lambda[sum(nnz<=optionsList$selfselect)])[-1]
        removeVar <- 
          c(removeVar,(possibleVar[-removeVar])[which(beta==0)])
      }
    }
    if(length(removeVar)>0) 
      possibleVar <- possibleVar[-removeVar]
    
    
    parents <- numeric(0)
    if(length(possibleVar)>1){
      gl <- glmnet::cv.glmnet(X[ ,possibleVar,drop=FALSE], 
                              as.numeric( X[,parentsOf[k]]),
                              intercept=TRUE)
      beta <- as.numeric(coef(gl))[-1]
      parents <- possibleVar[which(beta!=0)]
    }else{
      beta <- 0
    }
    result[[k]] <- parents
    attr(result[[k]],"parentsOf") <- parentsOf[k]
    
    if(pointEst) 
      attr(result[[k]],"coefficients") <- beta[beta!=0]
    
  }
  
  if(directed){
    for (k in 1:length(parentsOf)){
      parentsVar <- result[[k]]
      
      if(confBound)
        coefsParentsVar <-  as.numeric(attr(result[[k]],"coefficients"))
      
      childVar <- as.numeric(attr(result[[k]],"parentsOf"))
      parentsVarOrig <- parentsVar
      for(p in 1:length(parentsVarOrig)){
        idxParentsOfParent <- which(sapply(result, 
                                           function(j) 
                                             is.element(parentsVarOrig[p], 
                                                        as.numeric(attr(j, "parentsOf")))))
        
        if(length(idxParentsOfParent) == 0)
          next
        
        parentsOfParentP <- result[[idxParentsOfParent]]
        
        if(confBound)
          coefsParentsOfParentP <- as.numeric(attr(result[[idxParentsOfParent]],"coefficients"))
        
        
        if(is.element(childVar, parentsOfParentP)){
          
          if(confBound){
            coefsParentsVar <-  coefsParentsVar[parentsVar != parentsVar[p]]
            coefsParentsOfParentP <- coefsParentsOfParentP[parentsOfParentP != childVar]
          }
          
          parentsVar <- parentsVar[parentsVar != parentsVar[p]]
          parentsOfParentP <- parentsOfParentP[parentsOfParentP != childVar]
          
          result[[k]] <- parentsVar
          attr(result[[k]],"parentsOf") <- parentsOf[k]
          
          result[[idxParentsOfParent]] <- parentsOfParentP
          attr(result[[idxParentsOfParent]],"parentsOf") <- parentsVarOrig[p]
          
          if(confBound){
            attr(result[[k]],"coefficients") <- coefsParentsVar
            attr(result[[idxParentsOfParent]],"coefficients") <- coefsParentsOfParentP
          }
          
        }
        
      }
    }
  }
  
  list(resList = result, resMat = NULL)
}