runPC <- function(X, parentsOf, alpha, variableSelMat, setOptions, directed, verbose, 
                   result, ...){
  
  dots <- list(...)
  if(length(dots) > 0){
    warning("options provided via '...' not taken")
  }
  
  # additional options for PC
  optionsList <- list("indepTest"=pcalg::gaussCItest, 
                      "labels"=as.character(1:ncol(X)),
                      "fixedEdges"=NULL,
                      "NAdelete"=TRUE, "m.max"=Inf, "u2pd" = "relaxed", 
                      "skel.method"= "stable", "conservative"=FALSE,
                      "maj.rule"=FALSE, "solve.confl"=FALSE,
                      "numCores"=1)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)

  if(is.null(optionsList$suffStat)){
    suffStat <- list(C = cor(X), n = nrow(X))
  }else{
    suffStat <- optionsList$suffStat
  }
  
  pc.fit <- pcalg::pc(suffStat = suffStat, 
                      indepTest = optionsList$indepTest, p=ncol(X), 
               alpha = alpha, 
               labels=optionsList$labels,
               fixedGaps= if(is.null(variableSelMat)) NULL else (!variableSelMat), 
               fixedEdges = optionsList$fixedEdges, 
               NAdelete= optionsList$NAdelete, m.max= optionsList$m.max, 
               u2pd=optionsList$u2pd, 
               skel.method= optionsList$skel.method, 
               conservative= optionsList$conservative, 
               maj.rule= optionsList$maj.rule, 
               solve.confl = optionsList$solve.confl, 
               numCores=optionsList$numCores,
               verbose= verbose)
  pcmat <- as(pc.fit@graph, "matrix") 
  if(directed){
    warning("Removing undirected edges from estimated adjacency matrix.")
    pcmat <- pcmat * (t(pcmat)==0)
  }
  
  result <- vector("list", length = length(parentsOf))
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- which(as.logical(pcmat[, parentsOf[k]]))
    attr(result[[k]],"parentsOf") <- parentsOf[k]
  }
  
  if(length(parentsOf) < ncol(X)){
    pcmat <- pcmat[,parentsOf]
  }
  
  list(resList = result, resMat = pcmat)
}