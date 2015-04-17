


getParentsStable <- function(X, environment, interventions= NULL, EV=1, nodewise=TRUE,threshold=0.75, sampleSettings=1/sqrt(2), sampleObservations=1/sqrt(2),  parentsOf=1:ncol(X), method= c("ICP","hiddenICP","hiddenICE","pc","lingam","ges","gies","cam","rfci")[1],  alpha=0.1, variableSelMat=NULL,  excludeTargetInterventions= TRUE, onlyObservationalData= FALSE,   setOptions = list(), warnings=TRUE, nsim=100 ){

    if(0){
         EV=1
         nodewise=TRUE
         threshold=0.75
         parentsOf=1:ncol(X)
         sampleSettings=1/sqrt(2)
         sampleObservations=1/sqrt(2)
             
        environment= NULL

         interventions= NULL
        alpha=0.1
        variableSelMat=NULL
        excludeTargetInterventions= TRUE
        onlyObservationalData= FALSE
        setOptions = list()
        warnings=FALSE
        nsim=100
    }
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
        resmat <- hiddenICE(X,environment, alpha=EV, threshold=threshold, nsim=nsim,  sampleSettings=1/sqrt(2), sampleObservations=1/sqrt(2), nodewise=nodewise)$AhatAdjacency

    }else{
        for (sim in 1:nsim){
            useSettings <- sample( uniqueSettings, drawE(subs))
            ind <- which(  environment %in% useSettings)
            useSamples <- sort(sample(ind, round(length(ind)*sampleObservations)))
            
            res <- getParents(X[useSamples,], parentsOf=parentsOf,   interventions= interventions[useSamples], environment= environment[useSamples], method= method,  alpha= alpha, variableSelMat=variableSelMat,  excludeTargetInterventions= excludeTargetInterventions, onlyObservationalData= onlyObservationalData,  returnAsList=FALSE, confBound=TRUE, setOptions = setOptions, warnings=warnings)
            diag(res) <- 0
            reskeep <- 0* as(res,"matrix")
            quse <- drawE( q)
            if( nodewise){
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
        rem <- which( resmat < threshold)
        if(length(rem)>0) resmat[rem] <- 0
        resmat <- round(100*resmat)
    }
     rownames(resmat) <- if(is.null(colnames(X))) as.character(1:ncol(X)) else colnames(X)
    colnames(resmat) <- rownames(resmat)[parentsOf]
    return(resmat)
   
}










