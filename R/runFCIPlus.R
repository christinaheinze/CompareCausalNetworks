runFCIPlus <- function(X, parentsOf, alpha, setOptions, 
                       directed, verbose, ...){
  
  dots <- list(...)
  if(length(dots) > 0){
    warning("options provided via '...' not taken")
  }
  
  # additional options for FCI
  optionsList <- list("indepTest"=pcalg::gaussCItest,
                      "alpha"=alpha)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  suffStat <- list(C = cor(X), n = nrow(X))
  fci.fit <- pcalg::fciPlus(suffStat, 
                         indepTest = optionsList$indepTest, 
                         alpha = alpha,
                         p=ncol(X), 
                         verbose= verbose )
  fcimat <- fci.fit@amat
  
  if(directed){ 
    stop("directed currently not implemented for fciplus.")
    warning("Removing undirected edges from estimated connectivity matrix.")
    
    # fcimat <- fcimat * (t(fcimat)==0) #TODO: fix
  }
  
  result <- vector("list", length = length(parentsOf))
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- which(as.logical(fcimat[, parentsOf[k]]))
    attr(result[[k]],"parentsOf") <- parentsOf[k]
  }
  
  if(length(parentsOf) < ncol(X)){
    fcimat <- fcimat[,parentsOf]
  }
  
  list(resList = result, resMat = fcimat)
}