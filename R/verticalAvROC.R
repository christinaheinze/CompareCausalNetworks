verticalAvROC <- function(nSamplesToDraw, ROCs, filterBy = NULL){
  
  idxRm <- which(sapply(ROCs, is.null))
  if(length(idxRm) > 0) ROCs <- ROCs[-idxRm] 
  
  if(!is.null(filterBy)){
    uniqueComb <- unique(sapply(filterBy, 
                                function(f) 
                                  lapply(ROCs, function(roc) unique(roc[,f]))))
    
    dfTmp <- NULL
    for(comb in 1:nrow(uniqueComb)){
      settings <- uniqueComb[comb,,drop=FALSE]
      rownames(settings) <- NULL
      
      # filter ROCs according to comb
      ROCsFiltIdx <- which(sapply(ROCs, 
                                  function(l) 
                                    all(l[,filterBy,drop=FALSE] == settings)))
      ROCsFilt <- ROCs[ROCsFiltIdx]
      nROCCurves <- length(ROCsFilt)
      # nEmpty <- sum(sapply(ROCsFilt, is.null))
      out <- numeric(0)
      
      fprSeq <- seq(0, 1, length = nSamplesToDraw)
      
      for(s in fprSeq){
        tprsum <- 0
        for(c in 1:nROCCurves){
          # if(is.null(ROCsFilt[[c]])) next
          tprsum <- tprsum + tprForFpr(s, ROCsFilt[[c]])
        }
        out <- c(out, tprsum/nROCCurves)
      }
      # dfToAdd <- data.frame(fpr = fprSeq, tpr = out, filterBy)
      
      dfTmp <- rbind(dfTmp, data.frame(fpr = fprSeq, tpr = out, settings))
    }
    
    toReturn <- dfTmp
  }else{
    nROCCurves <- length(ROCs)
    # nEmpty <- sum(sapply(ROCs, is.null))
    out <- numeric(0)
    
    fprSeq <- seq(0, 1, length = nSamplesToDraw)
    
    for(s in fprSeq){
      tprsum <- 0
      for(c in 1:nROCCurves){
        # if(is.null(ROCs[[c]])) next
        tprsum <- tprsum + tprForFpr(s, ROCs[[c]])
      }
      out <- c(out, tprsum/nROCCurves)
    }
    toReturn <- data.frame(fpr = fprSeq, tpr = out)
  }
  toReturn
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

ROCdfAllMethods <- function(evalList, queries, nSamplesToDraw, filterBy = NULL){
  methods <- unique(unlist(lapply(evalList, function(l) names(l$evaluation) )))
  queriesList <- vector("list", length = length(queries))
  names(queriesList) <- queries
  
  for(q in queries){
    dfTmp <- NULL
    for(m in methods){
      if(is.null(filterBy)){
        runList <- lapply(evalList, function(l) l$evaluation[[m]][[q]])
        
        fprSeq <- seq(0,1,length = nSamplesToDraw)
        
        TPRAv <- verticalAvROC(nSamplesToDraw, runList)
        dfTmpM <- data.frame(method = m, TPRAv)
        
        dfTmp <- rbind(dfTmp, dfTmpM)
      }else{
        # for(filt in 1:length(filterBy)){
          
          runList <- lapply(evalList, function(l){
            tmp <- l$evaluation[[m]][[q]]
            if(!is.null(tmp)){
              givenConfig <- l$configs[filterBy]
              dfConf <- sapply(givenConfig, function(o) rep(o, times = nrow(tmp)))
              tmp <- data.frame(l$evaluation[[m]][[q]], dfConf)
            }
            tmp
          })
          
          TPRAv <- verticalAvROC(nSamplesToDraw, runList, filterBy)
          
          dfTmpM <- data.frame(method = m, TPRAv)
          
          dfTmp <- rbind(dfTmp, dfTmpM)
        # }
      }
      
    }
    queriesList[[q]] <- dfTmp
    
  }
  queriesList
}
