runFCI <- function(X, suffStat, parentsOf, alpha, variableSelMat, setOptions, directed, verbose, 
                    ...){
  
  dots <- list(...)
  if(length(dots) > 0){
    warning("options provided via '...' not taken")
  }
  
  # additional options for FCI
  optionsList <- list("indepTest"=pcalg::gaussCItest, 
                      "labels"=as.character(1:ncol(X)),
                      "skel.method"="stable", 
                      "type"="normal",
                      "fixedEdges"=NULL,
                      "NAdelete"=TRUE, 
                      "m.max"=Inf,
                      "pdsep.max"=Inf,
                      "rules"=rep(TRUE,10),
                      "doPdsep"=TRUE,
                      "biCC"=FALSE,
                      "conservative"=FALSE, 
                      "maj.rule"=FALSE, 
                      "numCores"=1)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  if(is.null(suffStat)){
    suffStat <- list(C = cor(X), n = nrow(X))
  }
  fci.fit <- pcalg::fci(suffStat, 
                         indepTest = optionsList$indepTest, 
                         alpha = alpha,
                         labels=optionsList$labels,
                         p=ncol(X), 
                         skel.method= optionsList$skel.method, 
                         type=optionsList$type,
                         fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), 
                         fixedEdges=optionsList$fixedEdges, 
                         NAdelete=optionsList$NAdelete, 
                         m.max=optionsList$m.max, 
                         pdsep.max=optionsList$pdsep.max, 
                         rules=optionsList$rules,
                         doPdsep=optionsList$doPdsep,
                         biCC=optionsList$biCC,
                         conservative= optionsList$conservative, 
                         maj.rule=optionsList$maj.rule,
                         numCores=optionsList$numCores,
                         verbose=verbose)
  fcimat <- fci.fit@amat
  
  if(is.logical(fcimat)){
    fcimat[fcimat] <- 1
    fcimat[!fcimat] <- 0
  }
  
  if(directed){ 
    stop("directed currently not implemented for fci.")
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