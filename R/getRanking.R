#' Estimate the connectivity matrix of a causal graph using stability selection.
#'
#' @description Estimates the connectivity matrix of a directed causal graph, 
#'  using various possible methods. Supported methods at the moment are backShift, 
#'  bivariateANM, bivariateCAM, CAM, hiddenICP, ICP, GES, GIES, LINGAM, 
#'  PC, regression and RFCI.
#'  Uses stability selection to select an appropriate sparseness.
#'
#' @param X A (nxp)-data matrix with n observations of p variables.
#' @param environment An optional vector of length n, where the entry for 
#' observation i is an index for the environment in which observation i took 
#' place (simplest case entries \code{1} for observational data and entries
#'  \code{2} for interventional data of unspecified type). Is required for 
#'  methods \code{ICP}, \code{hiddenICP}, \code{backShift}.
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
#' @param EV A bound on the expected number of falsely selected edges.
#' @param nodewise If \code{FALSE}, stability selection retains for each 
#'  subsample the largest overall entries in the connectivity matrix. 
#'  If \code{TRUE}, values are ordered row- and node-wise first and then the 
#'  largest entries in each row and column are retained. Error control is 
#'  valid (under exchangeability assumption) in both cases. The latter setting 
#'  \code{TRUE} is perhaps more robust and is the default. 
#' @param threshold The empirical selection frequency in (0.5,1) under 
#'  subsampling that needs to be surpassed for an edge to be selected.
#' @param nsim The number of resamples for stability selection.
#' @param sampleSettings The fraction of different environments to resample 
#'  in each resampling (at least two different environments will be selected so 
#'  the argument is without effect if there are just two different environments 
#'  in total).
#' @param sampleObservations The fraction of samples to resample in each 
#'  environment.
#' @param parentsOf The variables for which we would like to estimate the 
#' parents. Default are all variables.
#' @param method A string that specfies the method to use. The methods 
#' \code{pc} (PC-algorithm), \code{LINGAM} (LINGAM), \code{arges} (Adaptively 
#' restricted greedy equivalence search), \code{ges} 
#' (Greedy equivalence search), \code{gies} (Greedy interventional equivalence 
#' search),  \code{fci} (Fast causal inference)  
#' and \code{rfci} (Really fast causal inference) are imported from the 
#' package "pcalg" and are documented there in more detail, including the 
#' additional options that can be supplied via \code{setOptions}. The method 
#' \code{CAM} (Causal additive models) is documented in the package "CAM" and 
#' the methods \code{ICP} (Invariant causal prediction), \code{hiddenICP} 
#' (Invariant causal prediction with hidden variables) are from the package 
#' "InvariantCausalPrediction". The method \code{backShift} comes from the 
#' package "backShift". The method \code{mmhc} comes from the 
#' package "bnlearn". 
#' Finally, the methods \code{bivariateANM} and 
#' \code{bivariateCAM} are for now implemented internally but will hopefully 
#' be part of another package at some point in the near future.
#' @param alpha The level at which tests are done. This leads to confidence 
#' intervals for \code{ICP} and \code{hiddenICP} and is used internally for 
#' \code{pc} and \code{rfci}.
#' @param variableSelMat An optional logical matrix of dimension (pxp). An 
#' entry \code{TRUE} for entry (i,j) says that variable i should be considered 
#' as a potential parent for variable j and vice versa for \code{FALSE}. If the 
#' default value of \code{NULL} is used, all variables will be considered, but 
#' this can be very slow, especially for methods \code{pc}, \code{ges}, 
#' \code{gies}, \code{rfci} and \code{CAM}.
#' @param excludeTargetInterventions When looking for parents of variable k 
#' in 1,...,p, set to \code{TRUE} if observations where an intervention on 
#' variable k occured should be excluded. Default is \code{TRUE}.
#' @param onlyObservationalData If set to \code{TRUE}, only observational data 
#' is used. It will take the index in \code{environment} specified by 
#' \code{indexObservationalData}. If \code{environment} is \code{NULL}, all 
#' observations are used. Default is \code{FALSE}.
#' @param indexObservationalData Index in \code{environment} that encodes 
#' observational data. Default is \code{1}.
#' @param setOptions A list that can take method-specific options; see the 
#' individual documentations of the methods for more options and their 
#' possible values.
#' @param verbose If \code{TRUE}, detailed output is provided.
#' 
#' @return A sparse matrix, where a 0 entry in (j,k) corresponds to an estimate 
#' of 'no edge' \code{j} -> \code{parentsOf[k]}. Entries between 0 and 100 
#' give the selection percentage of this edge over all resamples (set to 0 if
#' below critical threshold) and all non-zero values are considered as selected
#' edges.  
#' 
#' @references 
#' Stability selection (2010):  N. Meinshausen and P. Buhlmann, 
#' Journal of the Royal Statistical Society: Series B, 72, 417-473
#' 
#' @author Nicolai Meinshausen \email{meinshausen@@stat.math.ethz.ch}, Christina
#'  Heinze \email{heinze@@stat.math.ethz.ch}
#' 
#' @seealso \code{\link{getParents}} for the underlying point-estimate of 
#' the causal graph.
#'
#' @keywords Causality, Graph estimations
#'  
getRanking <- function(X, environment, interventions=NULL, 
                       queries = c("isParent", "isMaybeParent", "isNoParent",
                                   "isAncestor","isMaybeAncestor", "isNoAncestor"),
                       method= c("ICP", "hiddenICP", "backShift", "pc", 
                                 "LINGAM", "ges", "gies", "CAM", "fci", "rfci",
                                 "regression", "bivariateANM", 
                                 "bivariateCAM")[1],  
                       alpha=0.1, variableSelMat=NULL, 
                       excludeTargetInterventions=TRUE, 
                       onlyObservationalData=FALSE, 
                       indexObservationalData=NULL, 
                       setOptions=list(), 
                       nsim=100, 
                       sampleSettings=1/sqrt(2), 
                       sampleObservations=1/sqrt(2),
                       verbose=FALSE, ...){
  # number of variables
  p <- ncol(X)
  
  # find unique settings
  uniqueSettings <- unique(environment)
  
  # number of settings to draw in each round
  subs <- sampleSettings* length(uniqueSettings)
  
  # init result list
  simList <- vector("list", length = nsim)
  resList <- vector("list", length = length(queries))
  names(resList) <- queries
  resList <- lapply(resList, function(i) i <- matrix(0, ncol = p, nrow = p))
  for(i in 1:length(resList)) attr(resList[[i]], "name") <- queries[i]
  
  for (sim in 1:nsim){
    
    # TODO move sampling to new function
    
    # find observations to use in this round
    # only use observational data?
    if(onlyObservationalData){
      
      if(!is.null(indexObservationalData)){
        # if indexObservationalData is given, extract observational data
        ind <- which(environment %in% indexObservationalData)
        warning(paste("Will use only environment", 
                      indexObservationalData,
                      "(= observational data?) among the", 
                      length(uniqueSettings),
                      "given distinct environments for method", 
                      method, "(", length(ind),"observations)"))
      }else{
        # if indexObservationalData is set to NULL, use all data points
        ind <- 1:nrow(X)
        
        warning(paste("Will use all observations for method", method, 
                      "assuming all data points are observational data (", 
                      length(ind),"observations)"))
      }
      
      # sample from observational data
      useSamples <- sort(sample(ind, round(length(ind)*sampleObservations)))
    }else{
      # if onlyObservationalData is false, 
      # sample the settings to use and draw a balanced subsample
      useSettings <- sample(uniqueSettings, drawE(subs))
      useSamples <- NULL
      for(s in 1:length(useSettings)){
        ind <- which(  environment %in% useSettings[s])
        useSamples <- 
          c(useSamples, sort(sample(ind, round(length(ind)*sampleObservations))))
      }
    }
    
    # permutation of columns
    permuteCols <- sample(p)
    
    colnamesX <- if(is.null(colnames(X))) as.character(1:ncol(X)) else colnames(X)
    
    # run getParents with this subsample
    res <- try(getParents(X[useSamples,permuteCols], 
                      environment= environment[useSamples], 
                      interventions=interventions[useSamples],
                      parentsOf=1:p, 
                      method= method,  alpha= alpha, 
                      variableSelMat=variableSelMat,  
                      excludeTargetInterventions= excludeTargetInterventions, 
                      onlyObservationalData=onlyObservationalData, 
                      indexObservationalData=indexObservationalData, 
                      returnAsList=FALSE, pointConf=FALSE, 
                      setOptions = setOptions, verbose = verbose, ...))
    
    if(inherits(res, "try-error")){
      cat(paste("Error in method", method,
                ". Skipping subsampling iteration. getParents() returned the following error:", 
                geterrmessage()))
      next
    }
    
    # redo permutation of columns
    res <- res[order(permuteCols),order(permuteCols)]
    rownames(res) <- colnames(res) <- colnamesX

    simList[[sim]] <- res
    
    resultForQueries <- convertForRanking(res, queries, method = method)
    resList <- lapply(resList, function(r) r <- r + resultForQueries[attr(r, "name") == names(resultForQueries)][[1]])
   
  }
  
  resList <- lapply(resList, function(r) r/nsim*100)
  
  resList <- lapply(resList, function(resmat) {
    # add row and column names to result matrix
    rownames(resmat) <- colnames(resmat) <- colnamesX
    resmat
  })
  
  ranking <- lapply(resList, function(r){
    arrayInd(order(r, getVecTobreakTies(r, resList), decreasing = T), .dim = dim(r))
  })
  
  toReturn <- list(ranking = ranking,
                   resList = resList, 
                   simEstimates = simList)
  return(toReturn)
}

getVecTobreakTies <- function(m, histList){
  if(attr(m, "name") == "isParent"){
    vec <- histList$isMaybeParent
  }else if(attr(m, "name") == "isMaybeParent"){
    vec <- histList$isParent
  }else if(attr(m, "name") == "isNoParent"){
    vec <-  max(histList$isMaybeParent) - histList$isMaybeParent
  }else if(attr(m, "name") == "isAncestor"){
    vec <-histList$isMaybeAncestor
  }else if(attr(m, "name") == "isMaybeAncestor"){
    vec <- histList$isAncestor
  }else if(attr(m, "name") == "isNoAncestor"){
    vec <-  max(histList$isMaybeAncestor) - histList$isMaybeAncestor
  }else{
    stop("Query not supported.")
  }
  
  if(sum(vec) == 0) vec <- rnorm(ncol(m)^2)
  
  vec
}
