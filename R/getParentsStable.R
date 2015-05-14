


getParentsStable <- function(X, environment, interventions= NULL, EV=1, nodewise=TRUE,threshold=0.75, 
                             sampleSettings=1/sqrt(2), sampleObservations=1/sqrt(2),  parentsOf=1:ncol(X), 
                             method= c("ICP","hiddenICP","hiddenICE","pc","lingam","ges","gies","cam","rfci")[1],  
                             alpha=0.1, variableSelMat=NULL,  excludeTargetInterventions= TRUE, 
                             onlyObservationalData= FALSE, indexObservationalData = 1, setOptions = list(), warnings=TRUE, nsim=100 ){
    p <- ncol(X)
    resmat <- matrix(0,p,length(parentsOf))

    q <- if(!nodewise) sqrt(EV*(2*threshold-1)*(p^2-p)) else sqrt(EV*(2*threshold-1))
    if(q<1){
        qmin <- if(!nodewise) 1/ ((2*threshold-1)*(p^2-p)) else 1/ (2*threshold-1)
        stop(paste("Number of selected edges in each subsample is less than 1. Need to increase EV (expected number of false positives) to at least ", signif(qmin,3), " to have the chance of getting an interesting result", if(nodewise) " (or switch 'nodewise' to FALSE)." else "."),sep="")
    }
    if(sampleSettings<=0) stop("sampleSettings needs to be positive")
    if(sampleObservations<=0) stop("sampleObservations needs to be positive")
    if(sampleSettings>1) stop("sampleSettings needs to be at most 1")
    if(sampleObservations>1) stop("sampleObservations needs to be at most 1")
    
    uniqueSettings <- unique(environment)
    subs <- sampleSettings* length(uniqueSettings)
    drawE <- function(x){
        z <- floor(x)
        if( runif(1) <  x-z) z <- z+1
        z <- max(1,z)
        return(z)
    }
    if(method=="hiddenICE"){ ## 'hiddenICE' is doing internal subsampling already
      
      optionsList <- list("covariance"=TRUE, "threshold"=0.75, "nsim"=100,"sampleSettings"=1/sqrt(2),
                          "sampleObservations"=1/sqrt(2), "nodewise"=TRUE, "tolerance"=10^(-4), "baseSettingEnv" = 1)                
      availableOptions <- names(optionsList)
      changeOptions <- availableOptions[ availableOptions %in% names(setOptions)]
      if(length(changeOptions)>0){
        for (option in changeOptions) optionsList[[option]] <- setOptions[[option]]
      }
      resmat <- hiddenICE(X, environment, covariance=optionsList$covariance,  alpha=EV, threshold =threshold, 
                          nsim=nsim,sampleSettings=1/sqrt(2), sampleObservations=1/sqrt(2), nodewise=nodewise, 
                          tolerance=optionsList$tolerance, baseSettingEnv = optionsList$baseSettingEnv)$AhatAdjacency
    }else{
        for (sim in 1:nsim){
            if(onlyObservationalData){
              # if indexObservationalData is given, extract observational data and sample from these
              if(!is.null(indexObservationalData)){
                ind <- which(environment %in% indexObservationalData)
                useSamples <- sort(sample(ind, round(length(ind)*sampleObservations)))
                if(warnings) warning(paste("Will use only environment", indexObservationalData,"(= observational data?) among the", length(uniqueSettings),
                                           "given distinct environments for method", method, "(", length(ind),"observations)"))
              }
#               else{
#                 # if indexObservationalData is set to NULL, use all data points
#                 useSamples <- 1:nrow(X)
#                 if(warnings) warning(paste("Will use all observations for method", method, "assuming all data points are observational data
#                                            (", nrow(X),"observations)"))
#               }
            }else{
              # if onlyObservationalData is false, sample from settings
              useSettings <- sample( uniqueSettings, drawE(subs))
              ind <- which(  environment %in% useSettings)
              useSamples <- sort(sample(ind, round(length(ind)*sampleObservations)))
            }
            res <- getParents(X[useSamples,], parentsOf=parentsOf, interventions=interventions[useSamples], environment= environment[useSamples], method= method,  alpha= alpha, variableSelMat=variableSelMat,  excludeTargetInterventions= excludeTargetInterventions, onlyObservationalData= onlyObservationalData, indexObservationalData = indexObservationalData, returnAsList=FALSE, confBound=TRUE, setOptions = setOptions, warnings=warnings)
            diag(res) <- 0
            reskeep <- 0* as(res,"matrix")
            quse <- drawE( q)
            if(nodewise){
                for (k in 1:p){
                    wh <- order( abs( res[,k]), rnorm(p), decreasing=TRUE)[1:quse]
                    reskeep[wh,k] <- 1
                    wh <- order( abs( res[k,]), rnorm(p), decreasing=TRUE)[1:quse]
                    reskeep[k,wh] <- 1
                }
            }else{
                selected <- sort(abs(res), decreasing =  TRUE)[1:quse]
                indicesSelected <- which(abs(res) >= min(selected), arr.ind = TRUE)
                reskeep[indicesSelected] <- 1
            }
            
            resmat <- resmat + reskeep/nsim
        }
        rem <- which(resmat < threshold)
        if(length(rem)>0) resmat[rem] <- 0
        resmat <- round(100*resmat)
    }
    rownames(resmat) <- if(is.null(colnames(X))) as.character(1:ncol(X)) else colnames(X)
    colnames(resmat) <- rownames(resmat)[parentsOf]
    return(resmat)
}









