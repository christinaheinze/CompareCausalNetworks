checkMode <- function(mode, method){

  if(method %in% c("fci", "rankFci", "rfci", "fciplus")){
    if(grepl("Parent", mode)){
      emitWarningOrError <- TRUE
      msg <- paste("The method", method, "cannot be used to estimate 
                   parental relations (only ancestral relations).")
    }
  }
  
  # warning(msg)
}
 

getGraphMode <- function(method){
  if(method %in% c("fci", "rankFci", "rfci", "fciplus")){
    return("PAG")
  }else if(method %in% c("pc", "rankPc", "arges", "rankArges", "ges", 
                         "rankGes", "gies", "rankGies")){
    return("CPDAG")
  }else if(method %in%  c("directLINGAM", "LINGAM", "mmhc", "ICP", "hiddenICP", "backShift",
                          "bivariateANM", "bivariateCAM", "CAM", "regression")){
    return("DAG")
  }else if(method %in% c("backShift")){
    return("CG")
  }else{
    stop(paste("Method", method, "not supported in getGraphMode()"))
  }
  
}

changeMode <- function(mode, method, result, p){
  if(mode == "raw") return(result)
  
  if(is.null(result$resMat)){
    mat <- matrix(0, ncol=p, nrow=p)
    for (k in 1:p){
      mat[result$resList[[k]],k] <- 1
    }
    result$resMat <- mat
  }
  
  resMat <- convertForRanking(amat = result$resMat, queries = mode, method = method)[[1]]

  result <- list()
  
  for (k in 1:p){
    result[[k]] <- which(resMat[, k] == 1)
    attr(result[[k]],"parentsOf") <- k
  }
  
  list(resList = result, resMat = resMat)
  
}