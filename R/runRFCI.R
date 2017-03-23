runRFCI <- function(X, parentsOf, alpha, variableSelMat, setOptions, directed, verbose, 
                    ...){
  
  dots <- list(...)
  if(length(dots) > 0){
    warning("options provided via '...' not taken")
  }
  
  # additional options for RFCI
  optionsList <- list("indepTest"=pcalg::gaussCItest,
                      "skel.method"="stable", "fixedEdges"=NULL,
                      "NAdelete"=TRUE, "m.max"=Inf,"rules"=rep(TRUE,10),
                      "conservative"=FALSE, "maj.rule"=FALSE)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  suffStat <- list(C = cor(X), n = nrow(X))
  rfci.fit <- pcalg::rfci(suffStat, indepTest = optionsList$indepTest, 
                   p=ncol(X), alpha=alpha, 
                   fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), 
                   fixedEdges=optionsList$fixedEdges, 
                   NAdelete=optionsList$NAdelete, m.max=optionsList$m.max, 
                   skel.method= optionsList$skel.method, 
                   conservative= optionsList$conservative, 
                   maj.rule=optionsList$maj.rule, rules=optionsList$rules, 
                   verbose= verbose )
  rfcimat <- rfci.fit@amat
  
  if(directed){ 
    stop("directed currently not implemented for fci.")
    warning("Removing undirected edges from estimated connectivity matrix.")
    
    # fcimat <- fcimat * (t(fcimat)==0) #TODO: fix
  }
  
  result <- vector("list", length = length(parentsOf))
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- which(as.logical(rfcimat[, parentsOf[k]]))
    attr(result[[k]],"parentsOf") <- parentsOf[k]
  }

  list(resList = result, resMat = rfcimat)
}