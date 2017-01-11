runARGES <- function(X, parentsOf, variableSelMat, setOptions, directed, verbose, 
                   result, ...){
  
  # additional options for GES
  if(is.null(setOptions$adaptive)){
    setOptions$adaptive <- "vstructures" #ARGES-CIG
  }else{
    if(setOptions$adaptive == "none")
      setOptions$adaptive <- "vstructures" #ARGES-CIG
  }
  
  runGES(X, parentsOf, variableSelMat, setOptions, directed, verbose, 
         result, ...)
}