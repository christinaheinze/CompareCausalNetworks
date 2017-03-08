verticalAvROCInterpolated <- function(nSamplesToDraw, ROCs, filterBy = NULL){
  # remove empty ROC curves
  idxRm <- which(sapply(ROCs, is.null))
  if(length(idxRm) > 0) ROCs <- ROCs[-idxRm] 
  
  # if filter is applied
  if(!is.null(filterBy)){
    # get unique combination of filter
    uniqueComb <- unique(sapply(filterBy, 
                                function(f) 
                                  lapply(ROCs, function(roc) unique(roc[,f]))))
    # for each unique combination
    dfTmpList <- lapply(1:nrow(uniqueComb), function(comb){
      # filter ROCs list by chosen combination
      settings <- uniqueComb[comb,,drop=FALSE]
      rownames(settings) <- NULL
      
      # filter ROCs according to comb; then average these curves
      ROCsFiltIdx <- which(sapply(ROCs, 
                                  function(l) 
                                    all(l[,filterBy,drop=FALSE] == settings)))
      ROCsFilt <- ROCs[ROCsFiltIdx]
      nROCCurves <- length(ROCsFilt)
      
      if(nROCCurves > 0){
        runsWithNAs <- sum(sapply(ROCsFilt, function(r) any(is.na(r$TPR))))
        merged.data.frame <- as.data.frame(data.table::rbindlist(ROCsFilt))
        
        # now ROCs is a list of data.frames all containing nSteps rows; average same values of FPR
        tmp <- lapply(filterBy, function(f)  merged.data.frame[,is.element(colnames(merged.data.frame), f)])
        names(tmp) <- filterBy
        tmp$FPR <- merged.data.frame$FPR
        meanROC <- aggregate(merged.data.frame$TPR, 
                             by = tmp, 
                             FUN = mean, na.rm = TRUE)
        
        sdROC <- aggregate(merged.data.frame$TPR, 
                             by = tmp, 
                             FUN = sd, na.rm = TRUE)
        
        tmp <- meanROC[,is.element(colnames(meanROC), filterBy),drop = FALSE]
        colnames(tmp) <- filterBy
        toReturnTmp <- data.frame(fpr = meanROC$FPR, 
                               tpr = meanROC$x, 
                               sdTpr = sdROC$x,
                               tmp,
                               nRuns = length(ROCsFilt) - runsWithNAs)
      }
      toReturnTmp
    })
    toReturn <- as.data.frame(data.table::rbindlist(dfTmpList))
  }else{
    nROCCurves <- length(ROCs)
    runsWithNAs <- sum(sapply(ROCs, function(r) any(is.na(r$TPR))))
    merged.data.frame <- as.data.frame(data.table::rbindlist(ROCs))
    # now ROCs is a list of data.frames all containing nSteps rows; average same values of FPR
    meanROC <- aggregate(merged.data.frame$TPR, 
                         by = list(FPR = merged.data.frame$FPR), 
                         FUN = mean, na.rm = TRUE)
    
    sdROC <- aggregate(merged.data.frame$TPR, 
                         by = list(FPR = merged.data.frame$FPR), 
                         FUN = sd, na.rm = TRUE)
    
    toReturn <- data.frame(fpr = meanROC$FPR, 
                           tpr = meanROC$x, 
                           sdTpr = sdROC$x,
                           nRuns = length(ROCs) - runsWithNAs)
  }
  toReturn
}

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
      
      if(nROCCurves > 0){
        out <- sout <- numeric(0)
        fprSeq <- seq(0, 1, length = nSamplesToDraw)
        for(s in fprSeq){
          tprsum <- 0
          nROCCurvesCorr <- nROCCurves

          for(c in 1:nROCCurves){
            increment <- tprForFpr(s, ROCsFilt[[c]])
            if(is.na(increment)){
              nROCCurvesCorr <- nROCCurvesCorr-1
            }else{
              tprsum <- tprsum + increment
            }
            
          }
          out <- c(out, if(nROCCurvesCorr > 0) tprsum/nROCCurvesCorr else NULL)
          if(nROCCurvesCorr == 0) sout <- c(sout, s)
        }
        if(length(sout) > 0) fprSeq <- fprSeq[!is.element(fprSeq, sout)]
        dfTmp <- rbind(dfTmp, data.frame(fpr = fprSeq, tpr = out, settings))
      }
      
    }
    
    toReturn <- dfTmp
  }else{
    nROCCurves <- length(ROCs)
    out <- sout <- numeric(0)
    fprSeq <- seq(0, 1, length = nSamplesToDraw)
    
    for(s in fprSeq){
      tprsum <- 0
      nROCCurvesCorr <- nROCCurves
      for(c in 1:nROCCurves){
        increment <- tprForFpr(s, ROCs[[c]])
        if(is.na(increment)){
            nROCCurvesCorr <- nROCCurvesCorr-1
          }else{
            tprsum <- tprsum + increment
          }
      }

      out <- c(out, if(nROCCurvesCorr > 0) tprsum/nROCCurvesCorr else NULL)
      if(nROCCurvesCorr == 0) sout <- c(sout, s)
      
    }
    if(length(sout) > 0) fprSeq <- fprSeq[!is.element(fprSeq, sout)]
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



fprForTpr <- function(tprSamp, ROC){
  idx <- which(ROC$TPR >= tprSamp)
  if(length(idx) == nrow(ROC)){
    res <- ROC[1, "FPR"]
  }else if(length(idx) == 0){
    res <- ROC[nrow(ROC), "FPR"]
  }else if(ROC[min(idx),"TPR"] == tprSamp){
    res <- ROC[min(idx), "FPR"]
  }else{
    idxMax <- min(idx)
    idxMaxP1 <- idxMax - 1
    slope <- (ROC[idxMaxP1, "FPR"] - ROC[idxMax, "FPR"])/
      (ROC[idxMaxP1, "TPR"] - ROC[idxMax, "TPR"])
    res <- ROC[idxMax, "FPR"] + slope*(tprSamp - ROC[idxMax, "TPR"])
  }
  res
}

tprForFprVec <- function(fprSamp, ROC){
  sapply(fprSamp, function(s) tprForFpr(s, ROC))
}



cutoff <- function(roc){
  if(any(is.na(roc$TPR))){
    dfTRet <- data.frame("FPRcut" = NA, "TPRcut" = NA)
  }else{
    dfTRet <- roc[which.min(abs(roc$FPR - (1 - roc$TPR))),c("FPR", "TPR")]
    colnames(dfTRet) <- c("FPRcut", "TPRcut")
    
  }
  dfTRet
}

computeAUC <- function(roc){
  if(any(is.na(roc$TPR))){
    return(NA)
  }else{
    # auc <- try(integrate(tprForFprVec, 0,1,roc, rel.tol = 0.01)$value)
    # if(inherits(auc, "try-error")){
    #   auc <- try(integrate(tprForFprVec, 0,1,roc, rel.tol = 0.1)$value)
    # }
    # auc
    mean(roc$TPR)
  }
}

ROCdfAllMethods <- function(evalList, queries, nSamplesToDraw, 
                            methods = unique(unlist(lapply(evalList, function(l) names(l$evaluation) ))),
                            filterBy = NULL,
                            filterByMethodOptions = FALSE,
                            methodOptions = NULL,
                            alreadyInterpolated = FALSE){
  queriesList <- vector("list", length = length(queries))
  names(queriesList) <- queries
  
  for(q in queries){
    # dfTmp <- NULL
    # for(m in methods){
    dfTmp <- lapply(methods, function(m){  
      runList <- lapply(evalList, function(l){
        tmp <- l$evaluation[[m]]
        if(!is.null(tmp)){
          if(!is.null(filterBy)) givenConfig <- l$configs[filterBy]
          toRet <- lapply(tmp, function(t){
            
            if(is.null(nrow(t$ROCs[[q]]))){
              if(!is.null(t$error)) if(t$error){
                cat(paste("\nMethod:", m, "with error", t$errorMsg, "run", l$hashConfigs))
              }
              return(NULL)
            } 
            
            givenOptions <- t$options
            if(!is.null(filterBy)) dfConf1 <- sapply(givenConfig, function(o) rep(o, times = nrow(t$ROCs[[q]])))
            dfConf2 <- sapply(givenOptions, function(o){
              optRep <- try(rep(o, times = nrow(t$ROCs[[q]])), silent = TRUE)
              if(inherits(optRep, "try-error")){
                nrowEvalDf <- nrow(t$ROCs[[q]])
                optRep <- rep(paste(as.character(o), collapse = ""), 
                              times = nrowEvalDf)
              }
              optRep
            })
            
            df <- data.frame(t$ROCs[[q]], dfConf2)
            
            df <- if(!is.null(filterBy)) data.frame(df, dfConf1) else df
            
            factors <- sapply(df, is.factor)
            df[,factors] <- sapply(df[,factors], as.character)
            df
            
            
          })
        }else{
          toRet <- NULL
        }
        toRet
      })
      runList <- unlist(runList, recursive = FALSE)  
      
      if(is.null(runList)) return(NULL)
      
      uniqueOpts <- unique(unlist(lapply(evalList, 
                                         function(l) lapply( l$evaluation[[m]], 
                                                             function(t) names(t$options)))))
      if(any(!is.element(methodOptions, uniqueOpts))){
        stop(paste("No results for methodOption", 
             methodOptions[which(!is.element(methodOptions, uniqueOpts))], 
             "in results for method", m))
      }
      if(alreadyInterpolated){
        TPRAv <- verticalAvROCInterpolated(nSamplesToDraw, runList, 
                               filterBy = if(filterByMethodOptions) c(methodOptions, filterBy) else filterBy)
      }else{
        TPRAv <- verticalAvROC(nSamplesToDraw, runList, 
                               filterBy = if(filterByMethodOptions) c(methodOptions, filterBy) else filterBy)
      }
      
      
      data.frame(method = m, TPRAv)
      
      # dfTmp <- rbind(dfTmp, dfTmpM)
    })

    queriesList[[q]] <- as.data.frame(data.table::rbindlist(dfTmp))
    
  }
  queriesList
}
