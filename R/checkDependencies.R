#' Checks whether the method is implemented and whether all dependencies are
#' installed
#'
#' @param method Method to use.
#'
checkDependencies <- function(method, fct = missingDependenciesMessage){
  if(method %in% c("ICP", "hiddenICP")){
    fct("InvariantCausalPrediction", method)
  }else if(method %in% c("rankPc", "pc", "LINGAM", "rankArges",
                         "arges", "rankGes", "ges", "gies", "rankGies",
                         "rankFci", "fci", "rfci", "fciplus")){
    fct("pcalg", method)
  }else if(method %in% c("rankArges", "arges")){
    fct("huge", method)
    fct("flare", method)
  }else if(method == "CAM"){
    fct("CAM", method)
  }else if(method == "regression"){
    fct("glmnet", method)
  }else if(method == "bivariateCAM"){
    fct("mgcv", method)
  }else if(method == "bivariateANM"){
    fct("kernlab", method)
    fct("mgcv", method)
  }else if(method == "mmhc"){
    fct("bnlearn", method)
  }else if(method == "directLINGAM"){
    fct("R.matlab", method)
  }else if(method == "backShift"){
    fct("backShift", method)
  }else{
    stop(paste("Method", method, "not (yet?) implemented."))
  }
}


missingDependenciesMessage <- function(package, method){
  if(!requireNamespace(package, quietly = TRUE)){
    stop(paste("The package '", package, "' is needed for ", method," to 
work. Please install it.", sep=""),
         call. = FALSE)
  }
}

checkRequireNamespace <- function(package, method){
  if(requireNamespace(package, quietly = TRUE)){
    TRUE
  }else{
    FALSE
  }
}