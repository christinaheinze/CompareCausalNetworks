runGIES <- function(X, interventions, parentsOf, variableSelMat, setOptions, 
                    directed, verbose, ...){

    # check validity of input arguments for GIES
  if(is.null(interventions)) 
    stop("'interventions' cannot be 'NULL' for method 'gies'")
  
  # additional optionsList for GIES
  optionsList <- list("phase"=c("forward", "backward", "turning"), 
                      "maxDegree"=integer(0), 
                      "lambda" = 0.5*log(nrow(X)))
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)

  isn <- which(sapply(interventions, is.null))
  if(length(isn)>0) { for (k in isn) interventions[[k]] <- numeric(0) }
  targets <- unique(interventions)
  target.index <- match(interventions, targets)
  targets <- as.list(lapply(targets, as.integer))
  
  if(any(table(unlist(targets)) == length(targets)))
    stop(paste("One variable is intervened on in all settings. At least one\n",
         "setting needs to be entirely different for 'gies' to run."), 
         call. = FALSE)
  
  score <- new("GaussL0penIntScore", 
               data=X, 
               targets=targets,  
               target.index = target.index,
               lambda = optionsList$lambda)
  

  tmp <- pcalg::gies(score, 
                     fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), 
                     phase=optionsList$phase, 
                     maxDegree=optionsList$maxDegree,
                     verbose=verbose, ...)
   
  giesmat <- as(tmp$essgraph, "matrix")
  giesmat[giesmat] <- 1
  giesmat[!giesmat] <- 0
  
  if(directed){
    warning("Removing undirected edges from estimated adjacency matrix.")
    giesmat <- giesmat * (t(giesmat)==0)
  }
  
  result <- vector("list", length = length(parentsOf))
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- which(giesmat[, parentsOf[k]] == 1)
    attr(result[[k]],"parentsOf") <- parentsOf[k]
  }
  
  if(length(parentsOf) < ncol(X)){
    giesmat <- giesmat[,parentsOf]
  }
  
  list(resList = result, resMat = giesmat)
}