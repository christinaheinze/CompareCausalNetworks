evaluateRanking <- function(trueDAG, estimatedRanking, queries){
  
  if(is.null(estimatedRanking)){
    return(NULL)
  }
  
  p <- ncol(trueDAG)
  
  # adjacency matrix of true DAG
  trueDAGAdj <- trueDAG
  trueDAGAdj[trueDAGAdj != 0] <- 1
  trueDAGAdj
  
  res <- lapply(queries, function(q) evaluateQuery(q, trueDAGAdj, estimatedRanking))
  names(res) <- queries
  res
}

evaluateQuery <- function(query, trueDAGAdj, estimatedRanking){
  p <- ncol(trueDAGAdj)
  idx <- which(names(estimatedRanking) == query)
  
  if(length(idx) == 0){
    stop(paste("Results for query `", query, "' not contained in results."))
  }else{
    est <- estimatedRanking[[idx]]
  }
  
  if(query == "isMaybeParent") query <- "isParent"
  if(query == "isMaybeAncestor") query <- "isAncestor"
  
  groundTruthConverted <- answerQueries(ancestralAmat = dag2ancestral(trueDAGAdj), 
                                        parentalAmat = trueDAGAdj, 
                                        query, 
                                        "DAG")
  getROCvalsVec(1:(p^2), est, groundTruthConverted)
}


# false positive rate
getROCvals <- function(k, est, groundTruth){
  p <- ncol(groundTruth)
  estTrunc <- est[1:k,]
  nPos <- sum(groundTruth)
  nNeg <- p^2 - sum(groundTruth)
  truePos <- sum(groundTruth[estTrunc])
  falsePos <- k - truePos
  c(FPR = falsePos/nNeg, TPR = truePos/nPos)
}

getROCvalsVec <- function(vec, est, truth){
  m <- sapply(vec, function(c) getROCvals(c, est, truth))
  data.frame(cbind(vec, t(m)))
}