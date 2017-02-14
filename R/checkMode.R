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
 
 
# changeMode <- function(mode, method, result){
#   if(mode == "raw")
#     return(result)
#   
#   if(method %in% c("fci", "rfci", "fciplus")){
#     if(mode == "ancestral"){
#       # how to handle different edge types:  
#       # a tail on an edge means that this tail is present in all MAGs in the 
#       # Markov equivalence class
#       #  
#       # an arrowhead on an edge means that this arrowhead is present in all 
#       # MAGs in the Markov equivalence class
#       # 
#       # a o-edgemark means that there is a at least one MAG in the Markov 
#       # equivalence class where the edgemark is a tail, and at least one 
#       # where the edgemark is an arrowhead
#       # 
#       # o-o
#       # o-
#       # o->
#       # ->
#       # <-> (hidden variables)
#       # - (selection variables)
#       # directed: 
#       # undirected
#       # bidirected
#     }
#   }else if(method %in% c("pc", "LINGAM", "arges", "ges", "gies", 
#                          "ICP", "hiddenICP", "backShift")){
#     # convert to ancestral
#     # CPDAG to PAG
#   }
# 
# }


getGraphMode <- function(method){
  if(method %in% c("fci", "rfci", "fciplus")){
    return("PAG")
  }else if(method %in% c("pc", "arges", "ges", "gies")){
    return("CPDAG")
  }else if(method %in%  c("LINGAM", "mmhc", "ICP", "hiddenICP", "backShift")){
    return("DAG")
  }else if(method %in% c("backShift")){
    return("CG")
  }
  
}