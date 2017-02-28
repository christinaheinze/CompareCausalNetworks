checkMode <- function(mode, method){
  emitWarningOrError <- FALSE
  
  if(method %in% c("fci", "rfci", "fciplus")){
    if(mode == "parental"){
      emitWarningOrError <- TRUE
      fct <- stop
      msg <- paste("The method", method, "cannot be used to estimate 
                   parental relations (only ancestral relations).")
    }
  }
  
  if(emitWarningOrError){
    fct(msg)
  }
}
 

getGraphMode <- function(method){
  if(method %in% c("fci", "rankFci", "rfci", "fciplus")){
    return("PAG")
  }else if(method %in% c("pc", "rankPc", "arges", "rankArges", "ges", 
                         "rankGes", "gies", "rankGies")){
    return("CPDAG")
  }else if(method %in%  c("directLINGAM", "LINGAM", "mmhc", "ICP", "hiddenICP", "backShift")){
    return("DAG")
  }else if(method %in% c("backShift")){
    return("CG")
  }else{
    stop(paste("Method", method, "not supported in getGraphMode()"))
  }
  
}