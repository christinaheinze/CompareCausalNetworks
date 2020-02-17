runRFCI <- function(X, parentsOf, alpha, variableSelMat, setOptions, directed, verbose, 
                    ...){
  
  dots <- list(...)
  if(length(dots) > 0){
    warning("options provided via '...' not taken")
  }
  
  # additional options for RFCI
  optionsList <- list("indepTest"=pcalg::gaussCItest, 
                      "labels"=as.character(1:ncol(X)),
                      "skel.method"="stable", "fixedEdges"=NULL,
                      "NAdelete"=TRUE, "m.max"=Inf,"rules"=rep(TRUE,10),
                      "conservative"=FALSE, "maj.rule"=FALSE,
                      "numCores"=1)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  suffStat <- list(C = cor(X), n = nrow(X))
  rfci.fit <- pcalg::rfci(suffStat, indepTest = optionsList$indepTest, 
                          alpha=alpha, 
                          labels=optionsList$labels,
                          p=ncol(X), 
                          fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), 
                          fixedEdges=optionsList$fixedEdges, 
                          NAdelete=optionsList$NAdelete, 
                          m.max=optionsList$m.max, 
                          skel.method= optionsList$skel.method, 
                          conservative= optionsList$conservative, 
                          maj.rule=optionsList$maj.rule, 
                          rules=optionsList$rules, 
                          numCores=optionsList$numCores,
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

  if(length(parentsOf) < ncol(X)){
    rfcimat <- rfcimat[,parentsOf]
  }
  
  list(resList = result, resMat = rfcimat)
}