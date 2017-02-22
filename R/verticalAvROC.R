verticalAvROC <- function(nSamplesToDraw, ROCs){
  nROCCurves <- length(ROCs)
  # lensCurves <- sapply(ROCs, nrow)
  out <- numeric(nSamplesToDraw)
  for(s in 1:nSamplesToDraw){
    tprsum <- 0
    for(c in 1:nROCCurves){
      tprsum <- tprsum + tprForFpr(s/nSamplesToDraw, ROCs[[c]])
    }
    out[s] <- tprsum/nROCCurves
  }
  out
}


tprForFpr <- function(fprSamp, ROC){
  idx <- which(ROC$FPR <= fprSamp)
  if(length(idx) == 0){
    res <- ROC[1, "TPR"]
  }else if(ROC[max(idx),"FPR"] == fprSamp){
    res <- ROC[max(idx), "TPR"]
  }else{
    if(max(idx) != nrow(ROC)){  # interpolate
      idxMax <- max(idx)
      idxMaxP1 <- idxMax + 1
      slope <- (ROC[idxMaxP1, "TPR"] - ROC[idxMax, "TPR"])/
        (ROC[idxMaxP1, "FPR"] - ROC[idxMax, "FPR"])
      res <- ROC[idxMax, "TPR"] + slope*(fprSamp - ROC[idxMax, "FPR"])
    }else{
      res <- ROC[max(idx), "TPR"]
    }
  }
  res
}

ROCdfAllMethods <- function(evalList, queries, nSamplesToDraw){
  methods <- unlist(unique(lapply(evalList, function(l) names(l$evaluation) )))
  queriesList <- vector("list", length = length(queries))
  names(queriesList) <- queries
  
  for(q in queries){
    dfTmp <- NULL
    for(m in methods){
      runList <- lapply(evalList, function(l) l$evaluation[[m]][[q]])
      # dfTmpM <- NULL
      
      fprSeq <- seq(0,1,length = nSamplesToDraw)
      TPRAv <- verticalAvROC(nSamplesToDraw, runList)
      dfTmpM <- data.frame(fpr = fprSeq, tpr = TPRAv, method = m)
      
      # runIDs <- names(runList)
      # for(r in 1:length(runList)){
      #   dfTmpM <- rbind(dfTmpM, data.frame(runList[[r]], runID = runIDs[r]))
      # }
      # dfTmpM$method <- m
      dfTmp <- rbind(dfTmp, dfTmpM)
    }
    queriesList[[q]] <- dfTmp
    
  }
  queriesList
}
