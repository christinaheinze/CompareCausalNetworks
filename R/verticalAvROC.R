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
      # 
      # dfToAdd <- data.frame(fpr = fprSeq, tpr = out)
      # for(i in 1:ncol(settings)){
      #   dfToAdd <- cbind(dfToAdd, rep(settings[,i], nrow(dfToAdd)))
      # }
      
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

ROCdfAllMethods <- function(evalList, queries, nSamplesToDraw, 
                            methods = unique(unlist(lapply(evalList, function(l) names(l$evaluation) ))),
                            filterBy = NULL,
                            filterByMethodOptions = FALSE,
                            methodOptions = NULL){
  queriesList <- vector("list", length = length(queries))
  names(queriesList) <- queries
  
  for(q in queries){
    dfTmp <- NULL
    for(m in methods){
      # if(is.null(filterBy)){
      #   
      #   runList <- lapply(evalList, function(l){
      #     tmp <- l$evaluation[[m]]
      #     lapply(tmp, function(t){
      #      if(!is.null(t)){
      #         givenOptions <- t$options
      #         dfConf <- sapply(givenOptions, function(o) rep(o, times = nrow(t[[q]])))
      #         data.frame(t[[q]], dfConf)
      #       }
      #     })
      #   })
      #   runList <- unlist(runList, recursive = FALSE)
      #   
      #   uniqueOpts <- unique(unlist(lapply(evalList, 
      #                                      function(l) lapply( l$evaluation[[m]], 
      #                                                          function(t) names(t$options)))))
      #   
      #   
      #   TPRAv <- verticalAvROC(nSamplesToDraw, runList, 
      #                          filterBy = if(filterByOptions) uniqueOpts else NULL)
      #   dfTmpM <- data.frame(method = m, TPRAv)
      #   
      #   dfTmp <- rbind(dfTmp, dfTmpM)
      # }else{
      #  
        runList <- lapply(evalList, function(l){
          tmp <- l$evaluation[[m]]
          if(!is.null(tmp)){
            if(!is.null(filterBy)) givenConfig <- l$configs[filterBy]
            toRet <- lapply(tmp, function(t){
              givenOptions <- t$options
              if(!is.null(filterBy)) dfConf1 <- sapply(givenConfig, function(o) rep(o, times = nrow(t[[q]])))
              dfConf2 <- sapply(givenOptions, function(o){
                optRep <- try(rep(o, times = nrow(t[[q]])), silent = TRUE)
                if(inherits(optRep, "try-error")){
                  optRep <- rep(paste(as.character(o), collapse = ""), times = nrow(t[[q]]))
                }
                optRep
              })
              
              df <- data.frame(t[[q]], dfConf2)
              
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
        
        if(is.null(runList)) next
        
        uniqueOpts <- unique(unlist(lapply(evalList, 
                                           function(l) lapply( l$evaluation[[m]], 
                                                               function(t) names(t$options)))))
        if(any(!is.element(methodOptions, uniqueOpts))){
          stop(paste("No results for methodOption", 
               methodOptions[which(!is.element(methodOptions, uniqueOpts))], 
               "in results for method", m))
        }
        TPRAv <- verticalAvROC(nSamplesToDraw, runList, 
                               filterBy = if(filterByMethodOptions) c(methodOptions, filterBy) else filterBy)
        
        dfTmpM <- data.frame(method = m, TPRAv)
        
        dfTmp <- rbind(dfTmp, dfTmpM)
      # }
      
    }
    queriesList[[q]] <- dfTmp
    
  }
  queriesList
}
