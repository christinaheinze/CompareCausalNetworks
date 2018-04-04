#' Helper function to change default options 
#'
#' @param availableOptions Options provided by function.
#' @param optionsToSet Options to be changed to user-provided values.
adjustOptions <- function(availableOptions, optionsToSet){
  
  namesAvailableOptions <- names(availableOptions)
  namesOptionsToSet <- names(optionsToSet)
  changeOptions <- namesAvailableOptions[namesAvailableOptions %in% namesOptionsToSet]
  if(length(changeOptions)>0){
    for (option in changeOptions) 
      availableOptions[[option]] <- optionsToSet[[option]]
  }
  
  additionalOptions <- namesOptionsToSet[!is.element(namesOptionsToSet, namesAvailableOptions)]
  if(length(additionalOptions)>0){
    for (option in additionalOptions) 
      availableOptions[[option]] <- optionsToSet[[option]]
  }
  
  availableOptions
}