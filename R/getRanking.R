#' Estimate a ranking of edges for causal relations in the underlying graph structure
#' using stability ranking.
#'
#' @description Estimates a ranking of edges for a given query, e.g. for 
#' parental relations in the underlying causal graph structure, using 
#' various possible methods. 
#' 
#' Supported methods at the moment are ARGES,
#' backShift, bivariateANM, bivariateCAM, CAM, FCI, FCI+, GES, GIES, hiddenICP, 
#' ICP, LINGAM, MMHC, rankARGES, rankFci, rankGES, rankGIES, rankPC, 
#' regression, RFCI and PC.
#'
#' @param X A \eqn{(n x p)}-data matrix with \eqn{n} observations of \eqn{p} variables.
#' @param environment A vector of length \eqn{n}, where the entry for 
#' observation \eqn{i} is an index for the environment in which observation \eqn{i} took 
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
#' @param queries One (or more of) "isParent", "isMaybeParent", "isNoParent",
#' "isAncestor","isMaybeAncestor", "isNoAncestor"
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
#' @param assumeNoSelectionVars Set to \code{TRUE} is you want to assume the absence 
#' of selection variables.
#' @param nsim The number of resamples for stability selection.
#' @param sampleSettings The fraction of different environments to resample 
#'  in each resampling (at least two different environments will be selected so 
#'  the argument is without effect if there are just two different environments 
#'  in total).
#' @param sampleObservations The fraction of samples to resample in each 
#'  environment.
#' @param verbose If \code{TRUE}, detailed output is provided.
#' @param ... Parameters to be passed to underlying method's function.
#' 
#' @return A list with the following entries:
#' \itemize{
#' \item \code{ranking} A list of length \code{length(queries)}. For each query,
#' the corresponding list entry contains a matrix of dimension \eqn{(p x p) x 2} 
#' with the ranking of edges. E.g. the first row indicates that the edge from 
#' ranking$isParent[1,1] to ranking$isParent[1,2] is the most likely edge according
#' to the method under consideration. 
#' \item \code{resList}  A list of length \code{length(queries)}. For each query,
#' the corresponding list entry contains a matrix  of dimension \eqn{(p x p)} with the counts for 
#' \eqn{A_{i,j} = 1} across the \code{nsim} subsamples.
#' \item \code{simEstimates} A list of length \code{nsim} with the method's 
#' output for each of the \code{nsim} subsamples.
#' }
#'
#' @details For both parental and ancestral relations, three queries are supported. 
#' The existence of a relation is assessed by the queries \code{isParent} and 
#' \code{isAncestor}; the absence of a relation is assessed by the queries 
#' \code{isNoParent} and \code{isNoAncestor}; the potential existence of a 
#' relation is addressed by the queries \code{isMaybeParent} and 
#' \code{isMaybeAncestor}.
#' 
#' All queries return a connectivity matrix which we denote by \eqn{A}. 
#' The interpretation of the entries of \eqn{A} differs according to the considered query:
#' 
#' \strong{Parental relations:} Queries concerning parental relations can only 
#' be answered by those methods under consideration that return a DAG, a CPDAG 
#' or a directed cyclic graph. When we say that a particular method cannot 
#' answer a given query, then the method's output with respect to this query 
#' will be the zero matrix. However, the eventual ranking for such a query will 
#' not necessarily be random due to the tie breaking scheme that is applied 
#' when ranking pairs of variables (see below).
#' 
#' \enumerate{
#' \item \code{isParent} In the connectivity matrix \eqn{A} returned by this 
#' query, the entry \eqn{A_{i,j} = 1} means that there is \emph{a directed edge} 
#' from node \eqn{i} to node \eqn{j} in the graph structure estimated by the 
#' method under consideration. Otherwise, \eqn{A_{i,j} = 0}.
#' \item \code{isMaybeParent} \eqn{A_{i,j} = 1} means that there is 
#' \emph{a directed or an undirected edge} from node \eqn{i} to node \eqn{j}
#' in the estimated graph structure. Otherwise, \eqn{A_{i,j} = 0}.
#' \item \code{isNoParent} \eqn{A_{i,j} = 1} means that there is neither a 
#' directed nor an undirected edge from node \eqn{i} to node \eqn{j} in the 
#' estimated graph structure. Otherwise, \eqn{A_{i,j} = 0}.
#' }
#' 
#' \strong{Ancestral relations:} Queries concerning ancestral relations can be 
#' answered by all methods under consideration.
#' \enumerate{
#' \item \code{isAncestor} \eqn{A_{i,j} = 1} means that there is a 
#' \emph{directed path} from node \eqn{i} to node \eqn{j} in the estimated graph 
#' structure. Otherwise, \eqn{A_{i,j} = 0}. In case of PAGs, directed paths can 
#' contain the edge types \eqn{i --> j} and \eqn{i --o j}. Including the latter 
#' edge type in this category implies that we exclude the existence of selection 
#' variables.
#' \item \code{isMaybeAncestor} \eqn{A_{i,j} = 1} then means that there is a 
#' path from node \eqn{i} to node \eqn{j} that contains directed and/or undirected 
#' edges. Otherwise, \eqn{A_{i,j} = 0}.	For PAGs, such paths can contain the edge 
#' types \eqn{i --> j}, \eqn{i --o j}, \eqn{i o-o j} and/or 
#' \eqn{i o-> j}. Otherwise, \eqn{A_{i,j} = 0}.
#' \item \code{isNoAncestor} \eqn{A_{i,j} = 1} means that there is neither a 
#' directed path nor a partially directed path from node \eqn{i} to node \eqn{j}
#'  in the estimated graph structure. Otherwise, \eqn{A_{i,j} = 0}.
#'  }
#' 
#' \strong{Stability ranking:} To obtain a ranking of edges for a given set of 
#' queries, we run the method under consideration on \code{nsims} random 
#' subsamples of the data. In each round, we draw samples from a fraction of 
#' settings, where the size of the fraction is specified by \code{sampleSettings}. 
#' In each chosen setting, we sample a fraction of observations 
#' uniformly at random without replacement, where the size of the fraction is 
#' specified by  \code{sampleObservations}. 
#' 
#' For each subsample we randomly 
#' permute the order of the variables in the input. 
#' Methods that are order-dependent can therefore not exploit any potential 
#' advantage stemming from a data matrix with columns ordered according to the 
#' causal ordering or a similar one. We then run the method on each subsample. 
#' 
#' For each subsample and a particular query, we obtain the corresponding 
#' connectivity matrix \eqn{A}. We can then rank all pairs of nodes \eqn{i,j} 
#' according to the frequency of the occurrence of \eqn{A_{i,j} = 1} across 
#' subsamples. Ties between pairs of variables can be broken with the results 
#' of the other queries if they are also computed as specified by \code{queries}; 
#' otherwise ties are broken at random:
#' 
#' \itemize{
#' \item If the query is \code{isParent}, ties are broken with counts for 
#' \code{isMaybeParent}. 
#' \item For the query \code{isMaybeParent} ties are broken with counts for 
#' \code{isParent}, i.e. in case of equal counts we give a preference to the 
#' edge that was considered more often to be a 'certain' parent. For methods 
#' returning DAGs this scheme makes the ranking for \code{isMaybeParent} equal 
#' to the result for \code{isParent}, up to the random tie breaking that is 
#' applied for \code{isParent}.
#' \item If the query is \code{isNoParent}, ties are broken according to which 
#' edge was selected less often in the query \code{isMaybeParent}. 
#' \item If the query is \code{isAncestor}, ties are broken with counts for 
#' \code{isMaybeAncestor}. 
#' \item For the query \code{isMaybeAncestor} ties are broken with counts 
#' for \code{isAncestor}, i.e. in case of equal counts we give a preference 
#' to the edge that was considered more often to be a 'certain' ancestor. 
#' For methods returning DAGs this scheme makes the ranking for \code{isMaybeAncestor} 
#' equal to the result for \code{isAncestor}, up to the random tie breaking 
#' that is applied for \code{isAncestor}.
#' \item If the query is \code{isNoAncestor}, ties are broken according to 
#' which one was selected less often in the query \code{isMaybeAncestor}. 
#' } 
#' 
#' If the tie breaking matrix defined according to these rules is 0, 
#' a matrix with standard normal random entries is used to break ties. 
#' Similarly, if there are remaining ties after applying the tie breaking rules 
#' described above, ties are broken randomly.
#'     
#' @author Christina Heinze-Deml \email{heinzedeml@@stat.math.ethz.ch}
#' 
#' @seealso \code{\link{getParents}} for the underlying point-estimate of 
#' the causal graph.
#'
#' @examples 
#' data("simDataInv")
#' X <- simDataInv$X
#' set.seed(1)
#' if(require(pcalg)){
#'   rank <- getRanking(X,
#'                 environment = simDataInv$environment,
#'                 queries = c("isParent","isMaybeParent"),
#'                 method = c("LINGAM"),
#'                 verbose = FALSE)
#'   # estimated ranking
#'   print(rank$ranking$isParent)
#'  
#'   # true adjacency matrix
#'   print(simDataInv$configs$trueA)
#' }else{
#'   cat("\nThe packages 'pcalg' is needed for the example to
#' work. Please install it.")
#' }
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
                       assumeNoSelectionVars = TRUE,
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
                      method= method,  
                      alpha= alpha, 
                      mode = "raw",
                      variableSelMat=variableSelMat,  
                      excludeTargetInterventions= excludeTargetInterventions, 
                      onlyObservationalData=onlyObservationalData, 
                      indexObservationalData=indexObservationalData, 
                      returnAsList=FALSE, 
                      sparse = FALSE,
                      directed=FALSE, 
                      pointConf = FALSE, 
                      setOptions = setOptions, 
                      assumeNoSelectionVars = assumeNoSelectionVars,
                      verbose = verbose, ...))
    
    if(inherits(res, "try-error")){
      errorMsg <- geterrmessage()
      
      cat(paste("Error in method", method,
                            ". getParents() returned the following error:",
                errorMsg))
      
      # system is computationally singular
      if(length(grep("system is computationally singular", errorMsg)) > 0){
        res <- matrix(0, nrow = p, ncol = p)
      }else{
        stop(paste("Error in method", method,
                   ". getParents() returned the following error:",
                   errorMsg))
      }
      
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
  
  resListForRanking <- resList
  
  # if(is.element("isMaybeParent", queries) & is.element("isParent", queries)){
  #   resListForRanking$isMaybeParent <- resListForRanking$isMaybeParent + resListForRanking$isParent
  # }
  
  # if(is.element("isMaybeAncestor", queries) & is.element("isAncestor", queries)){
  #   resListForRanking$isMaybeAncestor <- resListForRanking$isMaybeAncestor + resListForRanking$isAncestor
  # }
  
  ranking <- lapply(resListForRanking, function(r){
    arrayInd(order(r, 
                   getVecTobreakTies(r, resListForRanking), # break ties with results for other queries
                   rnorm(ncol(r)^2), # if still ties, break randomly
                   decreasing = T), .dim = dim(r))
  })
  
  list(ranking = ranking,
       resList = resList, 
       simEstimates = simList)
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
  
  if(!is.null(vec)){
    if(sum(vec) == 0) vec <- rnorm(ncol(m)^2)
  }else{
    vec <- rnorm(ncol(m)^2)
  }
  
  vec
}
