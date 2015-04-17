
getParents <- function(X,  environment=NULL, parentsOf=1:ncol(X),interventions= NULL,method= c("hiddenICP","ICP","hiddenICE","pc","lingam","ges","gies","cam","rfci","regression","bivariateANM","bivariateCAM")[1],  alpha=0.1, variableSelMat=NULL,  excludeTargetInterventions= TRUE, onlyObservationalData= FALSE,  returnAsList=FALSE, confBound=TRUE, setOptions = list(), warnings=TRUE){

    methodsList <- c("ICP","hiddenICP","hiddenICE","pc","lingam","ges","gies","cam","rfci","regression","bivariateANM","bivariateCAM")
    if(!method %in% methodsList){
        stop(paste("Method", method,"not (yet?) implemented"))
    }

    if( is.data.frame(X)) X <- as.matrix(X)
    if(!is.matrix(X)) stop("'X' needs to be a matrix")
    if( !all(as.numeric(parentsOf) %in% (1:ncol(X)))) stop("'parentsOf' needs to be a subset of 1:ncol(X)")
    if(!is.list(interventions) & !is.null(interventions)) stop("'interventions' needs to be a list or NULL")
    if(length(interventions)!=nrow(X) & !is.null(interventions)) stop("'interventions' needs to have as many entries as there are rows in 'X' (or be 'NULL')")
    if( is.null(environment) &is.null(interventions) & method %in% c("hiddenICP","ICP","hiddenICE","gies") ) stop(paste("'environment' and 'interventions' cannot both be 'NULL' for method", method))
    if( is.null(environment) ){
        environment <- match( interventions, unique(interventions))
        if((lu <- length(unique(environment)))>50) warning(paste("'environment' was set to NULL and has been created via  \n '> environment <- match( interventions, unique(interventions))'\n but this results in",lu," different environments (unique intervention combinations);\n very likely better to define a smaller number of environments using subject knowledge about the experiment by grouping various intervention targets into a single environment") )
    }
    if( length(environment) != nrow(X) & !is.null(environment)) stop("'environment' needs to have the same length as there are rows in 'X' (or be 'NULL')")
    if(alpha<0) stop("alpha needs to be positive")
    if(!is.null(variableSelMat)){
        if(!is.logical(variableSelMat)) stop("'variableSelMat' needs to be a matrix with boolean entries")
        if(nrow(variableSelMat)!=ncol(X)) stop("'variableSelMat' needs to have as many rows as there are variables (columns of 'X')")
    }
    
    
    if(onlyObservationalData ){ ## use only observational data
      if(length(ue <- unique(environment))>1){
        if(!is.null(interventions)){
          lengthinterventions <- numeric(length(ue))
          for (uc in 1:length(ue)){
            sel <- which( environment==ue[uc])
            lengthinterventions[uc] <- mean(sapply(interventions[sel],length))
          }
          usecl <- ue[which.min( lengthinterventions)]
        }else{
          usecl <- 1
        }
        if(warnings) warning(paste("will use only environment", ue[usecl],"(= observational data?) among the", length(ue),"given distinct environments for method", method))
        sel <- which(environment== ue[usecl])
        X <- X[sel,]
        environment <- environment[sel]
        interventions <- interventions[sel]
      }
    }else{
      if(excludeTargetInterventions & !is.null(interventions)){
        removeObsTarget <- list()
        for (parentsOfC in 1:length(parentsOf)){
          removeObsTarget[[parentsOfC]] <- which( sapply(interventions, function(x,a) a %in% x, a=parentsOf[parentsOfC]))
        }
      }
      
    }
    
    ## gather estimated parents of each "parentsOf" node in an elemnt of a list
    result <- list()
    for (k in 1:length(parentsOf)){
        result[[k]] <- numeric(0)
        attr(result[[k]],"parentsOf") <- parentsOf[k]
    }
    p <- ncol(X)
    n <- nrow(X)

    ## if package 'pcalg' is not loaded, still allow execution of code that is not making use of 'pcalg'
    functry <- try(gaussCItest,silent=TRUE)
    if(class(functry)=="try-error"){
      gaussCItest <- function(x) stop("function gaussCItest not loaded from package 'pcalg'")
    }
    functry <- try(selGam,silent=TRUE)
    if(class(functry)=="try-error"){
      selGam <- function(x) stop("function selGam not loaded from package 'pcalg'")
    }
    functry <- try(selGamBoost,silent=TRUE)
    if(class(functry)=="try-error"){
      selGamBoost <- function(x) stop("function selGamBoost not loaded from package 'pcalg'")
    }

    
    
    optionsList <- list("ICP" = list("gof" = 0.1,"test"="approximate","selection"="lasso","maxNoVariables"=7,"maxNoVariablesSimult"=3,"maxNoObs"=200,"stopIfEmpty"=TRUE, "selfselect"=NULL),
                        "hiddenICP"   = list("mode"="asymptotic", "selfselect"=NULL),
                        "hiddenICE"   = list("covariance"=TRUE, "threshold"=0.75, "nsim"=100,"sampleSettings"=1/sqrt(2),"sampleObservations"=1/sqrt(2), "nodewise"=TRUE),
                        "regression"   = list("selfselect"=NULL),
                        "gies"      = list("turning"=TRUE,"maxDegree"=integer(0),"verbose"=FALSE),
                        "ges"       = list("turning"=TRUE,"maxDegree"=integer(0),"verbose"=FALSE),
                        "pc"        = list("alpha" = 0.05, "indepTest" =gaussCItest, "fixedEdges"=NULL,"NAdelete"=TRUE,"m.max"=Inf,"u2pd","skel.method"= "stable","conservative"=FALSE,"maj.rule"=FALSE,"solve.confl"=FALSE,"verbose"=FALSE),
                        "lingam"    = list("output"=FALSE),
                        "cam"       = list("scoreName"="SEMGAM", "numCores"=1,  "output"=FALSE,"variableSel"=FALSE, "variableSelMethod"=selGamBoost,  "pruning"=FALSE, "pruneMethod"=selGam),
                        "rfci"      = list("alpha" = 0.05, "indepTest" =gaussCItest,"skel.method"="stable","fiedEdges"=NULL,"NAdelete"=TRUE,"m.max"=Inf,"rules"=rep(TRUE,10),"conservative"=FALSE,"maj.rule"=FALSE,"verbose"=FALSE),
                        "bivariateCAM"      = list("silent"=TRUE),
                        "bivariateANM"      = list("silent"=TRUE))
                        
    availableOptions <- names(optionsList[[method]])
    changeOptions <- availableOptions[ availableOptions %in% names(setOptions)]
    if(length(changeOptions)>0){
        for (option in changeOptions) optionsList[[method]][[option]] <- setOptions[[option]]
    }

    options <- optionsList[[method]]

    switch(method,
           "ICP" = {
                if( all((1:ncol(X)) %in% parentsOf) & is.null(interventions)) warning("ICP requires that no interventions occured on the target variables. \n In the current function call (a) all variables are considered as target variables (parentsOf=1:ncol(X)) and (b) the interventions are equal to NULL (and can thus not be removed for each variables). \n The results are likely misleading. Either target just specific variables by specifying 'parentsOf' or add the list where interventons occured (using argument 'interventions'.")
                
               for (k in 1:length(parentsOf)){
                   if(round(k/100)==(k/100)) cat(" ",k)

                 allobs <- 1:nrow(X)
                 if(excludeTargetInterventions & !is.null(interventions)){
                     if( length(removeObsTarget[[k]])>0) allobs <- allobs[ -removeObsTarget[[k]]]
                 }
                 possibleVar <- (1:ncol(X))
                 removeVar <- parentsOf[k]
                 if(!is.null(variableSelMat) & is.null(options$selfselect)){
                     if(ncol(variableSelMat)==length(parentsOf)) selc <- k else selc <- which( (1:ncol(X)) == parentsOf[k])
                     removeVar <- unique(c(removeVar,which( !variableSelMat[,selc] )))
                 }else{
                     if(!is.null(options$selfselect)){
                         
                         gl <- glmnet(X[, possibleVar[-removeVar],drop=FALSE], as.numeric(X[, parentsOf[k]]))
                         nnz <- apply(coef(gl)!=0, 2,sum)
                         beta <- coef(gl, s= gl$lambda[sum(nnz<=options$selfselect)])[-1]
                         removeVar <- c(removeVar,(possibleVar[-removeVar])[which(beta==0)])
                     }
                 }
                 if(length(removeVar)>0) possibleVar <- possibleVar[ -removeVar]

                 res <- ICP(X[allobs,possibleVar,drop=FALSE], as.numeric(X[allobs,parentsOf[k]]), ExpInd= environment[allobs], showAcceptedSets=FALSE,showCompletion=FALSE,alpha=alpha, gof=options$gof, test = options$test, selection =options$selection, maxNoVariables = options$maxNoVariables, maxNoVariablesSimult = options$maxNoVariablesSimult, maxNoObs=options$maxNoObs, stopIfEmpty = options$stopIfEmpty)
                 parents <- possibleVar[which( res$maximinCoefficients !=0)]
                 result[[k]] <- parents
                 if(confBound)  attr(result[[k]],"coefficients") <- res$maximinCoefficients[ which( res$maximinCoefficients !=0)]
               }
           },
           
           "hiddenICP" = {
               if( all((1:ncol(X)) %in% parentsOf) & is.null(interventions)) warning("hiddenICP requires that no interventions occured on the target variables. \n In the current function call (a) all variables are considered as target variables (parentsOf=1:ncol(X)) and (b) the interventions are equal to NULL (and can thus not be removed for each variables). \n The results are likely misleading. Either target just specific variables by specifying 'parentsOf' or add the list where interventons occured (using argument 'interventions'.")
               
               for (k in 1:length(parentsOf)){
                   if(round(k/100)==(k/100)) cat(" ",k)
                   allobs <- 1:nrow(X)
                   if(excludeTargetInterventions & !is.null(interventions)){
                       if(length(removeObsTarget[[k]])>0) allobs <- allobs[ -removeObsTarget[[k]]]
                   }
                   
                   possibleVar <- (1:ncol(X))
                   removeVar <- parentsOf[k]
                   if(!is.null(variableSelMat) & is.null(options$selfselect)){
                       if(ncol(variableSelMat)==length(parentsOf)) selc <- k else selc <- which( (1:ncol(X)) == parentsOf[k])
                       removeVar <- unique(c(removeVar,which( !variableSelMat[,selc] )))
                   }else{
                       if(!is.null(options$selfselect)){
                           
                           gl <- glmnet(X[, possibleVar[-removeVar],drop=FALSE], as.numeric(X[, parentsOf[k]]))
                           nnz <- apply(coef(gl)!=0, 2,sum)
                           beta <- coef(gl, s= gl$lambda[sum(nnz<=options$selfselect)])[-1]
                           removeVar <- c(removeVar,(possibleVar[-removeVar])[which(beta==0)])
                       }
                   }
                   if(length(removeVar)>0) possibleVar <- possibleVar[ -removeVar]
                   
                   res <- hiddenICP(X[allobs,possibleVar,drop=FALSE], as.numeric(X[allobs,parentsOf[k]]), environment[allobs],alpha=alpha,mode=options$mode)
                   parents <- possibleVar[wh <- which( res$maximinCoefficients !=0)]
                   result[[k]] <- parents
                   if(confBound)  attr(result[[k]],"coefficients") <- res$maximinCoefficients[ wh ]
               }
            },
           "hiddenICE" = {
               if( nrow(X) < ncol(X)) stop( "hiddenICE not suitable if there are more variables than observations")
               if( !is.null(variableSelMat)) warning( "option 'variableSelMat' not implemented for 'hiddenICE' -- using all variables")

               res <- hiddenICE(X, environment, covariance=options$covariance, alpha=alpha, threshold =options$threshold, nsim=options$nsim, sampleSettings=options$sampleSettings, sampleObservations=options$sampleObservations, nodewise=options$nodewise )
               for (k in 1:length(parentsOf)){
                   result[[k]] <- (wh <- which(res$AhatAdjacency[, k]!=0))
                   if(confBound)  attr(result[[k]],"coefficients") <- res$Ahat[ wh,k ]
               }
            },
           "regression" = {
               for (k in 1:length(parentsOf)){
                   if(round(k/100)==(k/100)) cat(" ",k)

                   possibleVar <- (1:ncol(X))
                   
                   possibleVar <- (1:ncol(X))
                   removeVar <- parentsOf[k]
                   if(!is.null(variableSelMat) & is.null(options$selfselect)){
                       if(ncol(variableSelMat)==length(parentsOf)) selc <- k else selc <- which( (1:ncol(X)) == parentsOf[k])
                       removeVar <- unique(c(removeVar,which( !variableSelMat[,selc] )))
                   }else{
                       if(!is.null(options$selfselect)){
                           
                           gl <- glmnet(X[, possibleVar[-removeVar],drop=FALSE], as.numeric(X[, parentsOf[k]]))
                           nnz <- apply(coef(gl)!=0, 2,sum)
                           beta <- coef(gl, s= gl$lambda[sum(nnz<=options$selfselect)])[-1]
                           removeVar <- c(removeVar,(possibleVar[-removeVar])[which(beta==0)])
                       }
                   }
                   if(length(removeVar)>0) possibleVar <- possibleVar[ -removeVar]
                   
                   
                   parents <- numeric(0)
                   if( length(possibleVar)>1){
                       gl <- cv.glmnet( X[ ,possibleVar,drop=FALSE], as.numeric( X[,parentsOf[k]]),intercept=TRUE)
                   
                       beta <- as.numeric(coef(gl))[-1]
                       parents <- possibleVar[which(beta!=0)]
                   }else{
                       beta <- 0
                   }
                   result[[k]] <- parents
                   if(confBound)  attr(result[[k]],"coefficients") <- beta[ beta!=0]

               }
           },
           "gies" = {
               isn <- which(sapply(interventions, is.null))
               if(length(isn)>0) { for (k in isn) interventions[[k]] <- numeric(0)}
               targets <- unique(interventions)
               target.index <- match(interventions, targets)

               score <- new("GaussL0penIntScore", data=X, targets =targets,  target.index = target.index)
               vers <- unlist(packageVersion('pcalg'))[1:2]
               if( all(vers >= c(2,1))){
                   tmp <- gies( score,fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat) )
               }else{
                   tmp <- gies( p, as.list(targets), score,fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat) )
               }
               giesmat <- as( tmp$essgraph, "matrix")
               for (k in 1:length(parentsOf)){
                   result[[k]] <- which(giesmat[, parentsOf[k]])
               }
           },
           "ges" = {
               vers <- unlist(packageVersion('pcalg'))[1:2]
               score <- new("GaussL0penObsScore", X)
               vers <- unlist(packageVersion('pcalg'))[1:2]
               if( all(vers >= c(2,1))){
                   G <- ges( score , fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), turning = options$turning, maxDegree=options$maxDegree, verbose=options$verbose)
               }else{
                   G <- ges( p,score , fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), turning = options$turning, maxDegree=options$maxDegree, verbose=options$verbose)
               }
               gesmat <- as(G$essgraph, "matrix")
               for (k in 1:length(parentsOf)){
                   result[[k]] <- which(gesmat[, parentsOf[k]])
               }
           },
           "pc" = {
               suffStat <- list(C = cor(X), n = nrow(X))
               pc.fit <- pc(suffStat, indepTest = options$indepTest, p = ncol(X), alpha = options$alpha, fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), fixedEdges = options$fixedEdges, NAdelete= options$NAdelete, m.max= options$m.max, u2pd=options$u2pd, skel.method= options$skel.method, conservative= options$conservative, maj.rule= options$maj.rule, solve.confl = options$solve.confl, verbose= options$verbose )
               pcmat <- as(pc.fit@graph, "matrix")
               for (k in 1:length(parentsOf)){
                   result[[k]] <- which(as.logical(pcmat[, parentsOf[k]]))
               }
           },
           "rfci" = {
               suffStat <- list(C = cor(X), n = nrow(X))
               rfci.fit <- rfci(suffStat, indepTest = options$indepTest, p = ncol(X), alpha = options$alpha, fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), fixedEdges = options$fixedEdges, NAdelete= options$NAdelete, m.max= options$m.max, skel.method= options$skel.method, conservative= options$conservative, maj.rule= options$maj.rule, rules = options$rules, verbose= options$verbose )
               rfcimat <- as(rfci.fit@amat, "matrix")
               for (k in 1:length(parentsOf)){
                   result[[k]] <- which(as.logical(rfcimat[, parentsOf[k]]))
               }
           },
           "lingam" = {
               if(nrow(X)<=ncol(X)) stop("LINGAM not suitable for high-dimensional data; need nrow(X) > ncol(X)")
              res <- LINGAM(X, output=options$output)
              lingammat <- res$Adj 
              for (k in 1:length(parentsOf)){
                   result[[k]] <- (wh <- which(lingammat[, parentsOf[k]]))
                   if(confBound)  attr(result[[k]],"coefficients") <- t(res$B)[ wh,parentsOf[k] ]
               }
           },
           "cam" = {
               if(!is.null(interventions)){
                   intervMat <- matrix(FALSE,nrow=nrow(X),ncol=ncol(X))
                   for (i in 1:length(interventions)){
                       if(length(interventions[[i]])>0) intervMat[i, interventions[[i]]] <- TRUE
                   }
                   cammat <- as(CAM(X,intervData=TRUE,intervMat=intervMat,variableSelMat=variableSelMat, scoreName=options$scoreName, numCores=options$numCores, output= options$output, variableSel=options$variableSel, variableSelMethod= options$variableSelMethod, pruning = options$pruning, pruneMethod=options$pruneMethod)$Adj,"matrix")
               }else{
                   cammat <- as(CAM(X)$Adj,"matrix")
               }
               for (k in 1:length(parentsOf)){
                   result[[k]] <- which(cammat[, parentsOf[k]]>0)
               }
           },
           "bivariateCAM" = {
               bivcammat <- bivariateCAM(X, parentsOf =parentsOf , variableSelMat = variableSelMat, silent = options$silent)
               for (k in 1:length(parentsOf)){
                   result[[k]] <- (wh <- which(bivcammat$causalParents[, parentsOf[k]]>0))
                   if(confBound)  attr(result[[k]],"coefficients") <- bivcammat$scoreMat[ wh,parentsOf[k] ]
               }
           },
           "bivariateANM" = {
               bivanmmat <- bivariateANM(X, parentsOf =parentsOf , variableSelMat = variableSelMat, silent = options$silent)
               for (k in 1:length(parentsOf)){
                   result[[k]] <- (wh <- which(bivanmmat$causalParents[, parentsOf[k]]>0))
                   if(confBound)  attr(result[[k]],"coefficients") <- bivanmmat$scoreMat[ wh,parentsOf[k] ]
               }
           },
           {
               warning(paste("method ", method," not implemented"))
           }
           )
    if(returnAsList){
        out <- result
    }else{
        rowind <- unlist(result)
        colind <- numeric(length(rowind))
        x <- unlist(lapply(result, function(x) attr(x,"coefficients")))
        if(is.null(x)) x <- 1
        cc <- 0
        for (k in 1:length(result)){
            norep <- length(result[[k]])
            if(norep>0){
                colind[ cc+(1:norep)] <- rep(k,norep)
                cc <- cc+norep
            }
        }
        resmat <- sparseMatrix(i=rowind,j=colind,x=x,dims=c(p, length(parentsOf)))
        colnames(resmat) <- parentsOf

        out <- resmat
    }
    rownames(out) <- if(is.null(colnames(X))) as.character(1:ncol(X)) else colnames(X)
    colnames(out) <- rownames(out)[parentsOf]
    
    return(out)
   
}










