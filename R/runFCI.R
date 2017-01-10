runFCI <- function(X, parentsOf, alpha, variableSelMat, setOptions, directed, verbose, 
                    result){
  
  # additional options for FCI
  optionsList <- list("indepTest"=pcalg::gaussCItest,
                      "skel.method"="stable", "fixedEdges"=NULL,
                      "NAdelete"=TRUE, "m.max"=Inf,"rules"=rep(TRUE,10),
                      "conservative"=FALSE, "maj.rule"=FALSE)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  suffStat <- list(C = cor(X), n = nrow(X))
  fci.fit <- pcalg::fci(suffStat, 
                         indepTest = optionsList$indepTest, 
                         alpha = alpha,
                         p=ncol(X), 
                         skel.method= optionsList$skel.method, 
                         fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), 
                         fixedEdges=optionsList$fixedEdges, 
                         NAdelete=optionsList$NAdelete, 
                         m.max=optionsList$m.max, 
                         rules=optionsList$rules,
                         conservative= optionsList$conservative, 
                         maj.rule=optionsList$maj.rule,
                         verbose= verbose )
  fcimat <- as(fci.fit@amat, "matrix")
  if(directed) fcimat <- fcimat * (t(fcimat)==0) #TODO: fix
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- which(as.logical(fcimat[, parentsOf[k]]))
  }
  
  result
}