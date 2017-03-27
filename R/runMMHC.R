runMMHC <- function(X, parentsOf, alpha, variableSelMat,
                    setOptions, verbose, ...){
  
  # additional options for MMHC
  optionsList <- list(whitelist = NULL, blacklist = NULL, test = NULL, score = NULL, 
                      k = log(nrow(X))/2,
                      B = NULL, restart = 0, perturb = 1, max.iter = Inf,
                      optimized = TRUE, strict = FALSE, debug = verbose)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  if(!is.null(variableSelMat)){
    blacklist <- as.data.frame(which(!variableSelMat, arr.ind = TRUE))
    colnames(blacklist) <- c("from", "to")
  }
         
  X <- as.data.frame(X)
  colnames(X) <- paste("V", 1:ncol(X), sep = "")
  
  res <- bnlearn::mmhc(X,
                       whitelist = optionsList$whitelist, 
                       blacklist = optionsList$blacklist, 
                       test = optionsList$test, 
                       score = optionsList$score,
                       alpha = alpha, 
                       B = optionsList$B, 
                       restart = optionsList$restart, 
                       perturb = optionsList$perturb, 
                       max.iter = optionsList$max.iter,
                       optimized = optionsList$optimized, 
                       strict = optionsList$strict, 
                       debug = optionsList$debug,
                       k = optionsList$k,
                       ...
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