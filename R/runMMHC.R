runMMHC <- function(X, parentsOf, alpha, variableSelMat,
                    setOptions, verbose, ...){
  
  
  dots <- list(...)
  if(length(dots) > 0){
    warning("options provided via '...' not taken")
  }
  
  # additional options for MMHC
  optionsList <- list(whitelist = NULL, blacklist = NULL, 
                      restrict.args = list(), 
                      maximize.args = list(),
                      debug = verbose)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  if(!is.null(variableSelMat)){
    blacklist <- as.data.frame(which(!variableSelMat, arr.ind = TRUE))
    colnames(blacklist) <- c("from", "to")
  }
         
  X <- as.data.frame(X)
  colnames(X) <- paste("V", 1:ncol(X), sep = "")
  
  res <- bnlearn::mmhc(x = X,
                       whitelist = optionsList$whitelist, 
                       blacklist = optionsList$blacklist, 
                       restrict.args = optionsList$restrict.args, 
                       maximize.args = optionsList$maximize.args,
                       debug = optionsList$debug
                       )
  
  # transform res to adjacency matrix
  mmhcmat <- sapply(res$nodes, function(node){
    parentsOfNode <- node$parents
    parents <- as.numeric(substr(parentsOfNode, start = 2, stop = nchar(parentsOfNode)))
    matrixRow <- rep(0, times = ncol(X))
    matrixRow[parents] <- 1
    matrixRow
  })
  
  result <- vector("list", length = length(parentsOf))
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- (wh <- which(mmhcmat[, parentsOf[k]] == 1))
    attr(result[[k]],"parentsOf") <- parentsOf[k]
  }
  
  if(length(parentsOf) < ncol(X)){
    mmhcmat <- mmhcmat[,parentsOf]
  }
   
  list(resList = result, resMat = mmhcmat)
}