evaluateRanking <- function(trueDAG, estimatedRanking, queries){
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
  
  groundTruthConverted <- answerQueries(ancestralAmat = dag2ancestral(trueDAGAdj), 
                                        parentalAmat = trueDAGAdj, 
                                        query, 
                                        "DAG")
  getROCvalsVec(1:(p^2), est, trueDAGAdj) #TODO 1:(p^2) or 1:(p^2)-p
}


# false positive rate
getROCvals <- function(k, est, trueDAGAdj, excludeDiag = FALSE){
  p <- ncol(trueDAGAdj)
  estTrunc <- est[1:k,]
  nPos <- sum(trueDAGAdj)
  nNeg <- (if(excludeDiag) (p^2-p) else p^2) - sum(trueDAGAdj)
  truePos <- sum(trueDAGAdj[estTrunc])
  falsePos <- k - truePos
  c(FPR = falsePos/nNeg, TPR = truePos/nPos)
}

getROCvalsVec <- function(vec, est, trueDAGAdj){
  m <- sapply(vec, function(c) getROCvals(c, est, trueDAGAdj))
  data.frame(cbind(vec, t(m)))
}