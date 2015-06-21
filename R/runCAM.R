runCAM<- function(X, interventions, parentsOf, variableSelMat, setOptions, 
                  directed, verbose, result){
  
  # additional options for CAM
  optionsList <- list("scoreName"="SEMGAM", "numCores"=1, "output"=FALSE, 
                      "variableSel"=FALSE, "variableSelMethod"=selGamBoost, 
                      "pruning"=FALSE, "pruneMethod"=selGam)
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
  if(!is.null(interventions)){
    intervMat <- matrix(FALSE,nrow=nrow(X),ncol=ncol(X))
    for (i in 1:length(interventions)){
      if(length(interventions[[i]])>0) intervMat[i, interventions[[i]]] <- TRUE
    }
    
    cammat <- as(CAM(X,intervData=TRUE,intervMat=intervMat,
                     scoreName=optionsList$scoreName, 
                     numCores=optionsList$numCores, 
                     output= optionsList$output, 
                     variableSel=optionsList$variableSel, 
                     variableSelMethod= optionsList$variableSelMethod, 
                     pruning = optionsList$pruning, 
                     pruneMethod=optionsList$pruneMethod)$Adj,"matrix")
  }else{
    cammat <- as(CAM(X)$Adj,"matrix")
  }
  if(directed) cammat <- cammat * (t(cammat)==0)
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- which(cammat[, parentsOf[k]]>0)
  }
  
  result
}