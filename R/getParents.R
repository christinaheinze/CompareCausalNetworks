#' Estimate the connectivity matrix of a causal graph
#'
#' @description Estimates the connectivity matrix of a directed causal graph, 
#' using various possible methods. Supported methods at the moment are PC, 
#' Lingam, GES, GIES, RFCI, CAM, and invariant prediction with or without hidden 
#' variables (ICP and hiddenICP) as well as backShift.
#'
#' @param X A (nxp)-data matrix with n observations of p variables.
#' @param environment An optional vector of length n, where the entry for 
#' observation i is an index for the environment in which observation i took 
#' place (simplest case entries \code{1} for observational data and entries
#'  \code{2} for interventional data of unspecified type). Is required for 
#'  methods "ICP", "hiddenICP", "backShift".
#' @param interventions A optional list of length n. The entry for observation
#'  i is a numeric vector that specifies the variables on which interventions 
#'  happened for observation i (a scalar if an intervention happened on just 
#'  one variable and \code{numeric(0)} if no intervention occured for this 
#'  observation). Is used for method \code{gies} but will generate the vector 
#'  \code{environment} if this is set to \code{NULL} (even though it might 
#'  generate too many different environments for some data so a hand-picked 
#'  vector \code{environment} is preferable). Is also used for \code{ICP} and 
#'  \code{hiddenICP} to exclude interventions on the target variable of 
#'  interest.
#' @param parentsOf The variables for which we would like to estimate the 
#' parents. Default are all variables.
#' @param method A string that specfies the method to use. The methods 
#' \code{pc} (PC-algorithm), \code{lingam} (Lingam), \code{ges} 
#' (Greedy equivalence search), \code{gies} (Greedy interventional equivalence 
#' search) and \code{rfci} (Really fast causal inference) are imported from the 
#' package "pcalg" and are documented there in more detail, including the 
#' additional options that can be supplied via \code{setOptions}. The method 
#' \code{cam} (Causal additive models) is documented in the package "cam" and 
#' the methods \code{ICP} (Invariant causal prediction), \code{hiddenICP} 
#' (Invariant causal prediction with hidden variables) are from the package 
#' "InvariantCausalPrediction".  The method \code{backShift} comes from the 
#' package "backShift". Finally, the methods \code{bivariateANM} and 
#' \code{bivariateCAM} are for now implemented internally but will hopefully 
#' be part of another package at some point in the near future.
#' @param alpha The level at which tests are done. This leads to confidence 
#' intervals for \code{ICP} and \code{hiddenICP} and is used internally for 
#' \code{pc} and \code{rfci}.
#' @param variableSelMat An optional logical matrix of dimension (pxp). An 
#' entry \code{TRUE} for entry (i,j) says that variable i should be considered 
#' as a potential parent for variable j and vice versa for \code{FALSE}. If the 
#' default value of \code{NULL} is used, all variables will be considered, but 
#' this can be very slow, especially for methods "pc", "ges", "gies", "rfci" 
#' and "cam".
#' @param excludeTargetInterventions When looking for parents of variable k 
#' in 1,...,p, set to \code{TRUE} if observations where an intervention on 
#' variable k occured should be excluded. Default is \code{TRUE}.
#' @param onlyObservationalData If set to \code{TRUE}, only observational data 
#' is used. It will take the index in \code{environment} specified by 
#' \code{indexObservationalData}. If \code{environment} is \code{NULL}, all 
#' observations are used. Default is \code{FALSE}.
#' @param indexObservationalData Index in \code{environment} that encodes 
#' observational data. Default is \code{1}.
#' @param returnAsList If set to \code{TRUE}, will return a list, where entry 
#' k is a list containing the estimated parents of variable k. The option 
#' \code{directed} will be ignored if set to \code{TRUE}. Default is \code{FALSE}.
#' @param confBound If \code{TRUE}, numerical estimates will be returned if 
#' possible. For methods \code{ICP} and \code{hiddenICP}, these are the values 
#' in the individual confidence intervals (at chosen level \code{alpha}) that 
#' are closest to 0; for other methods these are point estimates. Some methods 
#' do not return numerical point estimates; for these the output will remain 
#' binary 0/1 (no-edge/edge). Default is \code{FALSE}.
#' @param setOptions A list that can take method-specific options; see the 
#' individual documentations of the methods for more options and their 
#' possible values.
#' @param warnings
#' @param directed If \code{TRUE}, an edge will be returned if and only if an 
#' edge has been detected to be directed (ie entry will be set to 0 for entry 
#' (j,k) if both j->k and k-> j are estimated). Ignored if not the whole graph 
#' is estimated or if \code{returnAsList} is \code{TRUE}.
#'
#' @return If option \code{returnAsList} is \code{FALSE}, a sparse matrix, 
#' where a 0 entry in position (j,k) corresponds to an estimate of "no edge" 
#' \code{j} -> \code{parentsOf[k]}, while an entry 1 corresponds to an 
#' estimated egde. If option \code{confBound} is \code{TRUE}, the 1 entries 
#' will be replaced by numerical values that are either point estimates of the 
#' causal coefficients or confidence bounds (see above). 
#' If option \code{returnAsList} is \code{TRUE}, a list will be returned. 
#' The k-th entry in the list is the numeric vector with the indices of the 
#' estimated parents of node \code{parentsOf[k]}. 
#' 
#' @author Nicolai Meinshausen <meinshausen@stat.math.ethz.ch>, Christina
#' Heinze <heinze@stat.math.ethz.ch>
#' 
#' @seealso \code{\link{getParentsStable}} for stability selection-based 
#' estimation of the causal graph.
#' 
#' @examples
#' # 1st example:
#' # Simulate data with connectivity matrix A with assumptions 
#' # 1) hidden variables present
#' # 2) precise location of interventions is assumed unknown
#' # 3) different environments can be distinguished
# 
# simulate data
# set.seed(1)
# ## sample size n
# n <- 1000
# ## p=3 predictor variables and connectivity matrix A
# p  <- 3
# A <- diag(p)*0
# A[1,2] <- 0.8
# A[2,3] <- 0.8
# A[3,1] <- -0.4  
# 
# ## divide data in 10 different environments
# G <- 10
# environment <- rep(1:G, each=ceiling(n/G))[1:n]
# X <- Perturb <-  matrix(0,nrow=n,ncol=p)
# ## Input of hidden variables into each variable
# gamma <- rnorm(p)
# Input <- outer(W <- rnorm(n),gamma,FUN="*")
# ## simulate noise perturbations in each environment
# for (i in unique(environment)){
#   ind <- which(environment==i)
#   multiplier <- rexp(p)*3
#   Perturb[ind,] <- sweep(matrix(rnorm(length(ind)*p),ncol=p),
#                          2, multiplier,FUN="*")
# }
# ## iterate model to get stable solution
# ## (necessary only if feedbacks are included)
# niter <- 100
# for (iter in 1:niter){
#   X <- X \%*\% A + Input + Perturb 
# }
# 
# ####### apply all  methods given in vector 'methods'
# ####### (using all data pooled for pc/lingam/rfci --
# #######  -- can be changed with option 'onlyObservationalData=TRUE')
# methods <- c("hiddenICE", "lingam", "pc", "rfci","regression")
# 
# ## arrange graphical output into a rectangular grid
# sq <- ceiling(sqrt(length(methods)+1))
# par(mfrow=c(ceiling((length(methods)+1)/sq),sq))
# 
# ## plot and print true graph
# cat("\n true graph is  ------  \n" )
# print(A)
# plotGraph(A,main="true graph")
# 
# ## loop over all methods and compute and print/plot estimate
# for (method in methods){
#   cat("\n result for method", method,"  ------  \n" )
#   
#   ## Option 1): use this estimator as a point estimate
#   #   Ahat <- getParents(X, environment, method=method, alpha=0.1)
#   
#   ## Option 2): use a stability selection based estimator
#   ## with expected number of false positives bounded by EV=2
#   Ahat <- getParentsStable(X, environment, EV=2, method=method ,alpha=0.1)
#   
#   ## print and plot estimate
#   print(Ahat)
#   plotGraph(Ahat,main=paste("estimate for method",method))
# }
# 
# 
# 
# 
# 
# ##########################################
# ######## 2nd example:
# ######## Simulate data with connectivity matrix A with assumptions
# ######## 1) No hidden variables
# ######## 2) Precise location of interventions is known
# ##########################################
# ####### simulate data
# set.seed(1)
# ## sample size n
# n <- 2000
# ## p=5 predictor variables
# p  <- 5
# A <- diag(p)*0
# A[1,2] <- 0.8
# A[2,3] <- -0.8
# A[3,4] <- 0.8
# A[3,5] <- 0.8
# A[4,5] <- 0.3
# ## can add/remove feedback by using/not using
# A[5,2] <- 0.8 
# 
# ## choose explicity intervention targets
# interventions <- list()
# for (i in 1:n) interventions[[i]] <- sample(1:p,2)
# environment <- match(interventions, unique(interventions))
# X <- Perturb <-  matrix(0,nrow=n,ncol=p)
# ## Independent noise at each variable 
# Input <- matrix( rnorm(n*p),nrow=n)
# ## change level of noise for each intervention
# for (i in 1:n){
#   Perturb[i, interventions[[i]]] <- rnorm(length(interventions[[i]]))*5
# }
# ## iterate model to get stable solution 
# ## (necessary only if feedbacks are included)
# niter <- 100
# for (iter in 1:niter){
#   X <- X \%*\% A + Input + Perturb 
# }
# 
# 
# ####### apply possible  methods given in vector 'methods'
# ####### (using all data pooled for pc/lingam/rfci --
# #######    --can be changed with option 'onlyObservationalData=TRUE')
# methods <- c("hiddenICE","ICP","hiddenICP", "lingam", "pc", "rfci",
#              "regression","gies","ges")
# 
# ## arrange graphical output into a rectangular grid
# sq <- ceiling(sqrt(length(methods)+1))
# par(mfrow=c(ceiling((length(methods)+1)/sq),sq))
# 
# ## plot and print true graph
# cat("\n true graph is  ------  \n" )
# print(A)
# plotGraph(A,main="true graph")
# 
# ## loop over all methods and compute and print/plot estimate
# for (method in methods){
#   cat("\n result for method", method,"  ------  \n" )
#   
#   ## Option 1): use this estimator as a point estimate if desired:
#   #     Ahat <- getParents(X, environment, interventions=interventions,
#   #                       method=method ,alpha=0.1)
#   
#   ## Option 2): use a stability selection based estimator
#   ## with expected number of false positives bounded by EV=2
#   Ahat <- getParentsStable(X, environment,EV=2, interventions=interventions, 
#                            method=method ,alpha=0.1)
#   
#   
#   ## print and plot estimate
#   print(Ahat)
#   plotGraph(Ahat,main=paste("estimate for method",method))
# }
# 
# 
# 
# }
# #' 
# #' 
# #' 
# #' 
# 
# 

getParents <- function(X, environment = NULL, interventions = NULL, 
                       parentsOf = 1:ncol(X),
                       method = c("hiddenICP", "ICP", "backShift", "pc", 
                                  "lingam", "ges", "gies", "cam", "rfci", 
                                  "regression", "bivariateANM", 
                                  "bivariateCAM")[1],  
                       alpha = 0.1, variableSelMat = NULL,
                       excludeTargetInterventions = TRUE, 
                       onlyObservationalData = FALSE, indexObservationalData = 1,
                       returnAsList=FALSE, confBound = FALSE, 
                       setOptions = list(), warnings = TRUE, directed=TRUE){

    methodsList <- c("ICP","hiddenICP","backShift","pc","lingam","ges","gies","cam","rfci","regression","bivariateANM","bivariateCAM")
    if(!method %in% methodsList){
        stop(paste("Method", method,"not (yet?) implemented"))
    }

    if( is.data.frame(X)) X <- as.matrix(X)
    if(!is.matrix(X)) stop("'X' needs to be a matrix")
    if( !all(as.numeric(parentsOf) %in% (1:ncol(X)))) stop("'parentsOf' needs to be a subset of 1:ncol(X)")
    if(!is.list(interventions) & !is.null(interventions)) stop("'interventions' needs to be a list or NULL")
    if(length(interventions)!=nrow(X) & !is.null(interventions)) stop("'interventions' needs to have as many entries as there are rows in 'X' (or be 'NULL')")
    if( is.null(environment) &is.null(interventions) & method %in% c("hiddenICP","ICP","backShift","gies") ) stop(paste("'environment' and 'interventions' cannot both be 'NULL' for method", method))
    if( is.null(environment) & !is.null(interventions) ){
        environment <- match( interventions, unique(interventions))
        if((lu <- length(unique(environment)))>50) warning(paste("'environment' was set to NULL and has been created via  \n '> environment <- match( interventions, unique(interventions))'\n but this results in",lu," different environments (unique intervention combinations);\n very likely better to define a smaller number of environments using subject knowledge about the experiment by grouping various intervention targets into a single environment") )
    }
    if( length(environment) != nrow(X) & !is.null(environment)) stop("'environment' needs to have the same length as there are rows in 'X' (or be 'NULL')")
    if(alpha<0) stop("alpha needs to be positive")
    if(!is.null(variableSelMat)){
        if(!is.logical(variableSelMat)) stop("'variableSelMat' needs to be a matrix with boolean entries")
        if(nrow(variableSelMat)!=ncol(X)) stop("'variableSelMat' needs to have as many rows as there are variables (columns of 'X')")
    }
    if((length(unique(environment)) > 1) & onlyObservationalData & is.null(indexObservationalData)) stop("'indexObservationalData' needs to be specified")
    
    if(onlyObservationalData){ ## use only observational data
      if(length(ue <- unique(environment))>1){
#         if(!is.null(interventions)){
#           lengthinterventions <- numeric(length(ue))
#           for (uc in 1:length(ue)){
#             sel <- which( environment==ue[uc])
#             lengthinterventions[uc] <- mean(sapply(interventions[sel],length))
#           }
#           usecl <- ue[which.min( lengthinterventions)]
#         }else{
#           usecl <- indexObservationalData
#         }
#         if(warnings) warning(paste("will use only environment", ue[usecl],"(= observational data?) among the", length(ue),"given distinct environments for method", method))
#         sel <- which(environment== ue[usecl])
#         X <- X[sel,]
#         environment <- environment[sel]
#         interventions <- interventions[sel]
        
        sel <- which(environment== indexObservationalData)
        if(warnings) warning(paste("Will use only environment", indexObservationalData,"(= observational data?) among the", length(ue),
                                   "given distinct environments for method", method, "(", length(sel),"observations)"))
        
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
                        "backShift"   = list("covariance"=TRUE, "threshold"=0.75, "nsim"=100,"sampleSettings"=1/sqrt(2),"sampleObservations"=1/sqrt(2), "nodewise"=TRUE, "tolerance"=10^(-4), "baseSettingEnv" = 1),
                        "regression"   = list("selfselect"=NULL),
                        "gies"      = list("turning"=TRUE,"maxDegree"=integer(0),"verbose"=FALSE),
                        "ges"       = list("turning"=TRUE,"maxDegree"=integer(0),"verbose"=FALSE),
                        "pc"        = list("alpha" = 0.05, "indepTest" =gaussCItest, "fixedEdges"=NULL,"NAdelete"=TRUE,"m.max"=Inf,"u2pd","skel.method"= "stable","conservative"=FALSE,"maj.rule"=FALSE,"solve.confl"=FALSE,"verbose"=FALSE),
                        "lingam"    = list("verbose"=FALSE),
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
           "backShift" = {
               if( nrow(X) < ncol(X)) stop( "backShift not suitable if there are more variables than observations")
               if( !is.null(variableSelMat)) warning( "option 'variableSelMat' not implemented for 'backShift' -- using all variables")

               res <- try(backShift(X, environment, covariance=options$covariance, ev=alpha, 
                                threshold =options$threshold, nsim=options$nsim, sampleSettings=options$sampleSettings, 
                                sampleObservations=options$sampleObservations, nodewise=options$nodewise, tolerance=options$tolerance,
                                baseSettingEnv = options$baseSettingEnv), silent = FALSE)
               if(inherits(res, "try-error")){
                 cat("backShift -- no stable model found. Possible model mispecification. Returning the empty graph.\n")
                 res<- list(Ahat=0*diag(p), AhatAdjacency = 0*diag(p), varianceEnv = matrix(0, nrow = length(unique(environment)), ncol = p))
               }
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
               tryNewVersion <- try({ tmp <- gies( score,fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat) )},silent=TRUE)
               if(class(tryNewVersion)=="try-error"){
                   tmp <- gies( p, as.list(targets), score,fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat) )
               }
               giesmat <- as( tmp$essgraph, "matrix")
               if(directed) giesmat <- giesmat * (t(giesmat)==0)
               for (k in 1:length(parentsOf)){
                   result[[k]] <- which(giesmat[, parentsOf[k]])
               }
           },
           "ges" = {
               vers <- unlist(packageVersion('pcalg'))[1:2]
               score <- new("GaussL0penObsScore", X)
               vers <- unlist(packageVersion('pcalg'))[1:2]
               tryNewVersion <- try({G <- ges( score , fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), turning = options$turning, maxDegree=options$maxDegree, verbose=options$verbose)},silent=TRUE)
               if(class(tryNewVersion)=="try-error"){
                   G <- ges( p,score , fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), turning = options$turning, maxDegree=options$maxDegree, verbose=options$verbose)
               }
               gesmat <- as(G$essgraph, "matrix")
               if(directed) gesmat <- gesmat * (t(gesmat)==0)
               for (k in 1:length(parentsOf)){
                   result[[k]] <- which(gesmat[, parentsOf[k]])
               }
           },
           "pc" = {
               suffStat <- list(C = cor(X), n = nrow(X))
               pc.fit <- pc(suffStat, indepTest = options$indepTest, p = ncol(X), alpha = options$alpha, fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), fixedEdges = options$fixedEdges, NAdelete= options$NAdelete, m.max= options$m.max, u2pd=options$u2pd, skel.method= options$skel.method, conservative= options$conservative, maj.rule= options$maj.rule, solve.confl = options$solve.confl, verbose= options$verbose )
               pcmat <- as(pc.fit@graph, "matrix")
               if(directed) pcmat <- pcmat * (t(pcmat)==0)
  
               for (k in 1:length(parentsOf)){
                   result[[k]] <- which(as.logical(pcmat[, parentsOf[k]]))
               }
           },
           "rfci" = {
               suffStat <- list(C = cor(X), n = nrow(X))
               rfci.fit <- rfci(suffStat, indepTest = options$indepTest, p = ncol(X), alpha = options$alpha, fixedGaps=if(is.null(variableSelMat)) NULL else (!variableSelMat), fixedEdges = options$fixedEdges, NAdelete= options$NAdelete, m.max= options$m.max, skel.method= options$skel.method, conservative= options$conservative, maj.rule= options$maj.rule, rules = options$rules, verbose= options$verbose )
               rfcimat <- as(rfci.fit@amat, "matrix")
               if(directed) rfcimat <- rfcimat * (t(rfcimat)==0)

               for (k in 1:length(parentsOf)){
                   result[[k]] <- which(as.logical(rfcimat[, parentsOf[k]]))
               }
           },
           "lingam" = {
               if(nrow(X)<=ncol(X)) stop("LINGAM not suitable for high-dimensional data; need nrow(X) > ncol(X)")
              res <- LINGAM(X, verbose=options$verbose)
              lingammat <- res$Adj
               if(directed) lingammat <- lingammat * (t(lingammat)==0)

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
                                        #                    cammat <- as(CAM(X,intervData=TRUE,intervMat=intervMat,variableSelMat=variableSelMat, scoreName=options$scoreName, numCores=options$numCores, output= options$output, variableSel=options$variableSel, variableSelMethod= options$variableSelMethod, pruning = options$pruning, pruneMethod=options$pruneMethod)$Adj,"matrix")
                 cammat <- as(CAM(X,intervData=TRUE,intervMat=intervMat,scoreName=options$scoreName, numCores=options$numCores, output= options$output, variableSel=options$variableSel, variableSelMethod= options$variableSelMethod, pruning = options$pruning, pruneMethod=options$pruneMethod)$Adj,"matrix")
               }else{
                 cammat <- as(CAM(X)$Adj,"matrix")
               }
               if(directed) cammat <- cammat * (t(cammat)==0)

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










