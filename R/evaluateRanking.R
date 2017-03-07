evaluateRanking <- function(trueDAG, estimatedRanking, queries, 
                            interpolate = TRUE, nSteps = 100){
  
  if(is.null(estimatedRanking)){
    return(NULL)
  }
  
  p <- ncol(trueDAG)
  
  # adjacency matrix of true DAG
  trueDAGAdj <- trueDAG
  trueDAGAdj[trueDAGAdj != 0] <- 1
  trueDAGAdj
  
  res <- lapply(queries, function(q) evaluateQuery(q, trueDAGAdj, estimatedRanking, 
                                                   interpolate = interpolate, nSteps = nSteps))
  names(res) <- queries
  
  metrics <- lapply(queries, function(q){
    re <- res[[q]]
    cut <- cutoff(re)
    auc <- computeAUC(re)
    tprFpr0 <- tprForFpr(0, re)
    fprTpr1 <- fprForTpr(1, re)
    data.frame(cut, auc, tprFpr0, fprTpr1)
  }) 
  
  list(ROCs = res, metrics = metrics)
}

evaluateQuery <- function(query, trueDAGAdj, estimatedRanking, 
                          interpolate = TRUE, nSteps = 100){
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
  roc <- getROCvalsVec(0:(p^2), est, groundTruthConverted)
  if(interpolate){
    fprSeq <- seq(0,1,length.out = nSteps)
    tpr <- sapply(fprSeq, function(s) tprForFpr(s, roc))
    roc <- data.frame(vec = 0:(nSteps-1), FPR = fprSeq, TPR = tpr)         
  }
  roc
}


# false positive rate
getROCvals <- function(k, est, groundTruth){
  
  if(k == 0){
    return(c(FPR = 0, TPR = 0))
  }
  
  p <- ncol(groundTruth)
  estTrunc <- est[1:k,,drop = FALSE]
  nPos <- sum(groundTruth)
  nNeg <- p^2 - sum(groundTruth)
  truePos <- sum(groundTruth[estTrunc])
  falsePos <- k - truePos
  c(FPR = if(nNeg > 0) falsePos/nNeg else NA, 
    TPR = if(nPos > 0) truePos/nPos else NA)
}

getROCvalsVec <- function(vec, est, truth){
  m <- sapply(vec, function(c) getROCvals(c, est, truth))
  data.frame(cbind(vec, t(m)))
}