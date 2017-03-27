runDirectLINGAM <- function(X, parentsOf, pointEst, variableSelMat, setOptions, 
                            verbose, ...){
  
  dots <- list(...)
  matlab <- dots$matlab
  if(is.null(matlab)){
    stop("Need to provide matlab connection.")
  }
  
  if(is.character(matlab)){
    stop("Need to provide matlab connection. (is string)")
  }
  
  matlabFilesDir <- setOptions$matlabFilesDir
  if(is.null(matlabFilesDir)){
    stop("Need to provide directory with matlab files.")
  }
  
  if(nrow(X)<=ncol(X)) 
    stop("LINGAM not suitable for high-dimensional data; 
         need nrow(X) > ncol(X)")
 
  R.matlab::evaluate(matlab,paste("addpath(genpath('", matlabFilesDir,"'))", sep = ""))
  R.matlab::setVariable(matlab,data=t(X))
  R.matlab::evaluate(matlab,"Best=Dlingam(data);")
  lingammatCoef <- t(R.matlab::getVariable(matlab,"Best")$Best)
  lingammat <- lingammatCoef
  lingammat[lingammat != 0] <- 1
  
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