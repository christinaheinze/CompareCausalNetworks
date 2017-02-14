runGES <- function(X, parentsOf, variableSelMat, setOptions, directed, verbose, 
                   result, ...){
  
  # additional options for GES
  optionsList <- list("phases"= c("forward", "backward"),
                      "iterate"=FALSE,
                      "adaptive" = "none", 
                      "maxDegree"=integer(0))
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  score <- new("GaussL0penObsScore", X)
  G <- pcalg::ges(score, 
               fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), 
               adaptive = optionsList$adaptive,
               phase = optionsList$phase,
               iterate = optionsList$iterate,
               maxDegree=optionsList$maxDegree, 
               verbose=verbose, ...)
  gesmat <- as(G$essgraph, "matrix")
  gesmat[gesmat] <- 1
  gesmat[!gesmat] <- 0
  
  if(directed){
    warning("Removing undirected edges from estimated adjacency matrix.")
    gesmat <- gesmat * (t(gesmat)==0)
  }
  
  # for (k in 1:length(parentsOf)){
  #   result[[k]] <- which(gesmat[, parentsOf[k]] == 1) 
  # }
  
  gesmat
}