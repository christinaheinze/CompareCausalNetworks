runFCIPlus <- function(X, parentsOf, alpha, variableSelMat, setOptions, 
                       directed, verbose, result, ...){
  
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
  fcimat <- as(fci.fit@amat, "matrix")
  
  if(directed){ 
    stop("directed currently not implemented for fci.")
    # fcimat <- fcimat * (t(fcimat)==0) #TODO: fix
  }
  
  # for (k in 1:length(parentsOf)){
  #   result[[k]] <- which(as.logical(fcimat[, parentsOf[k]]))
  # }
  
  fcimat
}