runRESIT <- function(X, parentsOf, alpha, setOptions, 
                  verbose, ...){
  
  # additional options for CAM
  optionsList <- list("alpha" = alpha, "model" = train_gam, "parsModel" = list(), 
                      "indtest" = indtestHsic, "parsIndtest" = list(method = "ExactFastTrace"), 
                      "force_answer" = TRUE, "output" = verbose)
  
  dots <- list(...)
  if(length(dots) > 0){
    warning("options provided via '...' not taken")
  }
  
  # adjust according to setOptions if necessary
  optionsList <- adjustOptions(availableOptions = optionsList, 
                               optionsToSet = setOptions)
  
 adjmat <- RESIT(X, alpha = optionsList$alpha, model = optionsList$model, 
                 parsModel = optionsList$parsModel, 
                 indtest = optionsList$indtest, 
                 parsIndtest = optionsList$parsIndtest, 
                 force_answer = optionsList$force_answer, 
                 output = optionsList$output)

  result <- vector("list", length = length(parentsOf))
  
  for (k in 1:length(parentsOf)){
    result[[k]] <- which(adjmat[, parentsOf[k]]>0)
    attr(result[[k]],"parentsOf") <- parentsOf[k]
  }
  
  if(length(parentsOf) < ncol(X)){
    adjmat <- adjmat[,parentsOf]
  }
  
  list(resList = result, resMat = adjmat)
}