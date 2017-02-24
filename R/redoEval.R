redoEval <- function(simList, 
                     resList,
                     queries,
                     method,
                     redoConversion = FALSE){
  
  
  if(redoConversion){
    if(exists("resList")) rm(resList)
    resList <- vector("list", length = length(queries))
    names(resList) <- queries
    resList <- lapply(resList, function(i) i <- matrix(0, ncol = p, nrow = p))
    for(i in 1:length(resList)) attr(resList[[i]], "name") <- queries[i]
    
    nsim <- length(simList)
    
    
    for(s in 1:nsim){
      res <- simList[[s]]
      resultForQueries <- convertForRanking(res, queries, method = method)
      resList <- lapply(resList, function(r) r <- r + resultForQueries[attr(r, "name") == names(resultForQueries)][[1]])
    }
    
    resList <- lapply(resList, function(r) r/nsim*100)
    
    resList <- lapply(resList, function(resmat) {
      # add row and column names to result matrix
      rownames(resmat) <- colnames(resmat) <- colnamesX
      resmat
    })
    
  }

  ranking <- lapply(resList, function(r){
    arrayInd(order(r, getVecTobreakTies(r, resList), decreasing = T), .dim = dim(r))
  })

  list(ranking = ranking,
       resList = resList, 
       simEstimates = simList)
}