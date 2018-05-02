#' Estimate the connectivity matrix of a causal graph
#'
#' @description Estimates the connectivity matrix of a directed causal graph, 
#' using various possible methods. Supported methods at the moment are ARGES,
#' backShift, bivariateANM, bivariateCAM, CAM, FCI, FCI+, GES, GIES, hiddenICP, 
#' ICP, LINGAM, MMHC, rankARGES, rankFci, rankGES, rankGIES, rankPC, 
#' regression, RFCI and PC.
#' 
#' @param X A \eqn{(n} x \eqn{p)}-data matrix with n observations of  \eqn{p} variables.
#' @param environment An optional vector of length \eqn{n}, where the entry for 
#' observation \eqn{i} is an index for the environment in which observation \eqn{i} took 
#' place (Simplest case: entries \code{1} for observational data and entries
#'  \code{2} for interventional data of unspecified type. Encoding for observational
#'  data can be changed with \code{indexObservationalData}). Is required for 
#'  methods \code{ICP}, \code{hiddenICP} and \code{backShift}.
#' @param interventions A optional list of length \eqn{n}. The entry for observation
#'  \eqn{i} is a numeric vector that specifies the variables on which interventions 
#'  happened for observation \eqn{i} (a scalar if an intervention happened on just 
#'  one variable and \code{integer(0)} if no intervention occured for this 
#'  observation). Is used for methods \code{gies}, \code{rankGies} and \code{CAM} and will 
#'  generate the vector \code{environment} if the latter is set to \code{NULL}.
#'  (However, this might generate too many different environments for some data 
#'  sets, so a hand-picked vector \code{environment} is preferable). Is also used 
#'  for \code{ICP} and \code{hiddenICP} to exclude interventions on the target 
#'  variable of interest.
#' @param parentsOf The variables for which we would like to estimate the 
#' parents. Default are all variables. Currently only used with \code{mode = "raw"}. 
#' Speeds up computation for methods \code{bivariateANM},
#' \code{bivariateCAM}, \code{ICP}, \code{hiddenICP} and \code{regression}; 
#' for other methods only affects output. Also see \code{variableSelMat} for possibly
#' speeding up computational time by restricting the set of potential parents
#' for a variable. 
#' @param method A string that specfies the method to use. The methods 
#' \code{pc} (PC-algorithm), \code{LINGAM} (LINGAM), \code{arges} (Adaptively 
#' restricted greedy equivalence search), \code{ges} 
#' (Greedy equivalence search), \code{gies} (Greedy interventional equivalence 
#' search),  \code{fci} (Fast causal inference), \code{fciplus}  
#' and \code{rfci} (Really fast causal inference) are imported from the 
#' package "pcalg" and are documented there in more detail, including the 
#' additional options that can be supplied via \code{setOptions}. 
#' The "rank versions" of arges, fci, ges, gies and pc are based on [1]. The method 
#' \code{CAM} (Causal additive models) is documented in the package "CAM" and 
#' the methods \code{ICP} (Invariant causal prediction), \code{hiddenICP} 
#' (Invariant causal prediction with hidden variables) are from the package 
#' "InvariantCausalPrediction". The method \code{backShift} comes from the 
#' package "backShift". The method \code{mmhc} comes from the 
#' package "bnlearn". Finally, the methods \code{bivariateANM} and 
#' \code{bivariateCAM} are for now implemented internally but will hopefully 
#' be part of another package at some point in the near future.
#' @param alpha The level at which tests are done. This leads to confidence 
#' intervals for \code{ICP} and \code{hiddenICP} and is used internally for 
#' \code{pc}, \code{rankPc}, \code{mmhc}, \code{fci}, \code{rankFci}, \code{fciplus}
#'  and \code{rfci}. For all other methods \code{alpha} is not used.
#' @param mode Determines output type - can be "raw" or one of the queries "isParent", 
#' "isMaybeParent", "isNoParent", "isAncestor","isMaybeAncestor", "isNoAncestor".
#' If "raw", \code{getParents()} returns the connectivity matrix computed by the
#' specified method in sparse matrix format if \code{sparse} is set to \code{TRUE}; 
#' else in dense matrix format (or as list if \code{returnAsList = TRUE}). 
#' The options \code{directed} and  \code{pointConf} will be ignored for 
#' all modes except for "raw" if set to \code{TRUE}. The different mode types
#' are explained in the help for \code{\link{getRanking}}.
#' @param variableSelMat An optional logical matrix of dimension  \eqn{(p} x \eqn{p)}. An 
#' entry \code{TRUE} for entry \eqn{(i,j)} says that variable \eqn{i} should be considered 
#' as a potential parent for variable \eqn{j} and vice versa for \code{FALSE}. If the 
#' default value of \code{NULL} is used, all variables will be considered, but 
#' this can be very slow, especially for methods \code{pc}, \code{ges}, 
#' \code{gies}, \code{rfci} and \code{CAM}. Ignored for methods \code{backShift}, 
#' \code{fciplus}, \code{LINGAM} and \code{CAM}.
#' @param excludeTargetInterventions When looking for parents of variable \eqn{k} 
#' in \eqn{1,...,p}, set to \code{TRUE} if observations where an intervention on 
#' variable \eqn{k} occured should be excluded. Default is \code{TRUE}. Used
#' in  \code{ICP} and \code{hiddenICP}.
#' @param onlyObservationalData If set to \code{TRUE}, only observational data 
#' is used. It will take the index in \code{environment} specified by 
#' \code{indexObservationalData}. If \code{environment} is \code{NULL}, all 
#' observations are used. Default is \code{FALSE}.
#' @param indexObservationalData Index in \code{environment} that encodes 
#' observational data. Default is \code{1}.
#' @param returnAsList If set to \code{TRUE}, will return a list, where entry 
#' \eqn{k} is a list containing the estimated parents of variable \eqn{k}. 
#' Default is \code{FALSE}.
#' @param sparse If set to \code{TRUE} and \code{returnAsList} is \code{FALSE},
#' output matrix will be in sparse matrix format.
#' @param pointConf If \code{TRUE}, numerical estimates will be returned if 
#' possible. For methods \code{ICP} and \code{hiddenICP}, these are the values 
#' in the individual confidence intervals (at chosen level \code{alpha}) that 
#' are closest to 0; for other methods these are point estimates. Some methods 
#' do not return numerical point estimates; for these the output will remain 
#' binary 0/1 (no-edge/edge). Default is \code{FALSE}. Only supported in mode "raw".
#' @param setOptions A list that can take method-specific options; see the 
#' individual documentations of the methods for more options and their 
#' possible values.
#' @param directed If \code{TRUE}, an edge will be returned if and only if an 
#' edge has been detected to be directed. I.e. entry will be set to 0 for entry 
#' \eqn{(j,k)} if both \eqn{j -> k} and \eqn{k -> j} are estimated 
#' (\code{ICP}, \code{hiddenICP}, \code{regression}), if \eqn{j -- k} is undirected 
#' (in case of CPDAGs) or if the edge type is not of type \eqn{i --> j} in case of 
#' PAGs. If \code{assumeNoSelectionVars = TRUE} the edge type \eqn{i --o j} is also 
#' considered 'directed' for methods returning PAGs.  Default is \code{FALSE}. 
#' Only supported in mode "raw".
#' @param assumeNoSelectionVars Set to \code{TRUE} is you want to assume the absence 
#' of selection variables.
#' @param verbose If \code{TRUE}, detailed output is provided.
#' @param ... Parameters to be passed to underlying method's function.
#'
#' @return If option \code{returnAsList} is \code{FALSE}, a sparse matrix, 
#' where a 0 entry in position (j,k) corresponds to an estimate of "no edge" 
#' \code{j} -> \code{k}, while an entry 1 corresponds to an 
#' estimated egde. If option \code{pointConf} is \code{TRUE}, the 1 entries 
#' will be replaced by numerical values that are either point estimates of the 
#' causal coefficients or confidence bounds (see above). 
#' If option \code{returnAsList} is \code{TRUE}, a list will be returned. 
#' The k-th entry in the list is the numeric vector with the indices of the 
#' estimated parents of node \code{k}. 
#'  
#' @references 
#' \enumerate{
#' \item Naftali Harris and Mathias Drton: PC Algorithm for Nonparanormal 
#' Graphical Models. J. Mach. Learn. Res. 14(1) 2013.
#' }
#' 
#' @author Christina Heinze-Deml \email{heinzedeml@@stat.math.ethz.ch}, 
#'  Nicolai Meinshausen \email{meinshausen@@stat.math.ethz.ch}
#' 
#' @seealso \code{\link{getParentsStable}} for stability selection-based 
#' estimation of the causal graph.
#' 
#' @examples
#' ## load the backShift package for data generation and plotting functionality
#' if(require(backShift) & require(pcalg)){
#'   # Simulate data with connectivity matrix A with assumptions
#'   # 1) hidden variables present
#'   # 2) precise location of interventions is assumed unknown
#'   # 3) different environments can be distinguished
#'   
#'   ## simulate data
#'   myseed <- 1
#'   
#'   # sample size n
#'   n <- 10000
#'   
#'   # p=3 predictor variables and connectivity matrix A
#'   p <- 3
#'   labels <- c("1", "2", "3")
#'   A <- diag(p)*0
#'   A[1,2] <- 0.8
#'   A[2,3] <- 0.8
#'   A[3,1] <- -0.4
#'   
#'   # divide data in 10 different environments
#'   G <- 10
#'   
#'   # simulate
#'   simResult <- backShift::simulateInterventions(n, p, A, G, intervMultiplier = 3,
#'                noiseMult = 1, nonGauss = TRUE, hiddenVars = TRUE,
#'                knownInterventions = FALSE, fracVarInt = NULL, simulateObs = TRUE,
#'                seed = myseed)
#'   X <- simResult$X
#'   environment <- simResult$environment
#'   
#'   ## apply all  methods given in vector 'methods'
#'   ## (using all data pooled for pc/LINGAM/rfci/ges -- can be changed with option
#'   ## 'onlyObservationalData=TRUE')
#'   
#'   methods <- c("backShift", "LINGAM") #c("pc", "rfci", "ges")
#'   
#'   # select whether you want to run stability selection
#'   stability <- FALSE
#'   
#'   # arrange graphical output into a rectangular grid
#'   sq <- ceiling(sqrt(length(methods)+1))
#'   par(mfrow=c(ceiling((length(methods)+1)/sq),sq))
#'   
#'   ## plot and print true graph
#'   cat("\n true graph is  ------  \n" )
#'   print(A)
#'   plotGraphEdgeAttr(A, plotStabSelec = FALSE, labels = labels, thres.point = 0,
#'    main = "TRUE GRAPH")
#'   
#'   ## loop over all methods and compute and print/plot estimate
#'   for (method in methods){
#'     cat("\n result for method", method,"  ------  \n" )
#'   
#'     if(!stability){
#'       # Option 1): use this estimator as a point estimate
#'       Ahat <- getParents(X, environment, method=method, alpha=0.1, pointConf = TRUE)
#'     }else{
#'       # Option 2): use a stability selection based estimator
#'       # with expected number of false positives bounded by EV=2
#'       Ahat <- getParentsStable(X, environment, EV=2, method=method, alpha=0.1)
#'     }
#'   
#'    # print and plot estimate (point estimate thresholded if numerical estimates
#'    # are returned)
#'     print(Ahat)
#'     if(!stability)
#'       plotGraphEdgeAttr(Ahat, plotStabSelec = FALSE, labels = labels,
#'        thres.point = 0.05,
#'        main=paste("POINT ESTIMATE FOR METHOD\n", toupper(method)))
#'     else
#'       plotGraphEdgeAttr(Ahat, plotStabSelec = TRUE, labels = labels,
#'        thres.point = 0, main = paste("STABILITY SELECTION
#'        ESTIMATE\n FOR METHOD", toupper(method)))
#'    }
#' }else{
#'     cat("\nThe packages 'backShift' and 'pcalg' are needed for the examples to
#' work. Please install them.")
#' }
#'  
#' 
#'  
#' @keywords Causality, Graph estimations
#'  
getParents <- function(X, environment = NULL, interventions = NULL, 
                       parentsOf = 1:ncol(X),
                       method= c("arges", "backShift", "bivariateANM", 
                                 "bivariateCAM", "CAM", 
                                 "fci", "fciplus", "ges", "gies", "hiddenICP",
                                 "ICP", "LINGAM", "mmhc", "rankArges",  "rankFci",
                                 "rankGes", "rankGies", "rankPc", "rfci", "pc",
                                 "regression")[12],  
                       alpha = 0.1, 
                       mode = c("raw", "parental", "ancestral")[1],
                       variableSelMat = NULL,
                       excludeTargetInterventions = TRUE, 
                       onlyObservationalData = FALSE, 
                       indexObservationalData = 1,
                       returnAsList=FALSE, 
                       sparse = FALSE,
                       directed=FALSE, 
                       pointConf = FALSE, 
                       setOptions = list(), 
                       assumeNoSelectionVars = TRUE,
                       verbose = FALSE, ...){

    # check whether method is supported and dependencies are installed
    checkDependencies(method)
  
    # check whether mode is supported by method
    checkMode(mode, method)
    
    # check validity of other input arguments
    if(is.data.frame(X)) X <- as.matrix(X)
    if(!is.matrix(X)) stop("'X' needs to be a matrix")
    if(!all(as.numeric(parentsOf) %in% (1:ncol(X)))) 
      stop("'parentsOf' needs to be a subset of 1:ncol(X)")
    if(length(setdiff(1:ncol(X), parentsOf)) > 0 & mode != "raw") 
      stop("Combination of parentsOf and mode other than raw currently not implemented.")
    if(length(setdiff(1:ncol(X), parentsOf)) > 0 & !returnAsList) 
      stop("Combination of parentsOf and returnAsList=FALSE is currently not implemented.")
    if(!is.list(interventions) & !is.null(interventions)) 
      stop("'interventions' needs to be a list or NULL")
    if(length(interventions)!=nrow(X) & !is.null(interventions)) 
      stop("'interventions' needs to have as many entries as there are rows 
           in 'X' (or be 'NULL')")
    if(is.null(environment) &is.null(interventions) & method %in% 
       c("hiddenICP","ICP","backShift","gies") ) 
      stop(paste("'environment' and 'interventions' cannot 
                 both be 'NULL' for method", method))
    if(is.null(environment) & !is.null(interventions)){
        environment <- match(interventions, unique(interventions))
        if((lu <- length(unique(environment)))>50) 
          warning(paste(
                  "'environment' was set to NULL and has been created via \n 
                  '> environment <- match(interventions,unique(interventions))'\n 
                  but this results in", lu," different environments",
                  "(unique intervention combinations);\n 
                  very likely better to define a smaller number of environments 
                  using subject knowledge about the experiment by grouping 
                  various intervention targets into a single environment"))
    }
    if(length(environment) != nrow(X) & !is.null(environment)) 
      stop("'environment' needs to have the same length as there 
           are rows in 'X' (or be 'NULL')")
    if(!is.null(alpha)){ if(alpha<0) stop("alpha needs to be positive") }
    if(!is.null(variableSelMat)){
        if(!is.logical(variableSelMat)) 
          stop("'variableSelMat' needs to be a matrix with boolean entries")
        if(nrow(variableSelMat)!=ncol(X)) 
          stop("'variableSelMat' needs to have as many rows as there 
               are variables (columns of 'X')")
    }
    if(!is.null(variableSelMat) & method %in% c("backShift", "fciplus",
                                                "directLINGAM", "LINGAM", "CAM")) 
      warning(paste("option 'variableSelMat' not implemented for method",
            method, " -- using all variables"))
    if(directed & method %in% c("hiddenICP", "ICP", "regression", "LINGAM",
                                "mmhc", "CAM", "backShift")){
      warning("Option 'directed' is ignored for ICP, hiddenICP,
           LINGAM, mmhc, CAM, backShift and regression.")
    }
    
    if(directed & mode != "raw"){
      directed <- FALSE
      warning("Setting directed to FALSE. Only supported when mode is 'raw'.")
    }
    
    if(pointConf & mode != "raw"){
      pointConf <- FALSE
      warning("Setting pointConf to FALSE. Only supported when mode is 'raw'.")
    }
   
    # eval options
    setOptions <- lapply(setOptions, function(l) eval(l))
    
    # find unique settings
    uniqueSettings <- unique(environment)
    
    # select observations to use
    # use observational data only or 
    # use data from different environments/interventions
    if(onlyObservationalData & length(uniqueSettings) > 1){ 
      # use only observational data
      
      if(!is.null(indexObservationalData)){
        # if indexObservationalData is given, extract observational data
        sel <- which(environment %in% indexObservationalData)
        warning(paste("Will use only environment", 
                      indexObservationalData,
                      "(= observational data?) among the", 
                      length(uniqueSettings),
                      "given distinct environments for method", 
                      method, "(", length(sel),"observations)"))
      
      }else{
        # if indexObservationalData is set to NULL, use all data points
        sel <- 1:nrow(X)
        
        warning(paste("Will use all observations for method", method, 
                      "assuming all data points are observational data (", 
                      length(sel),"observations)"))
      }
      
      X <- X[sel,]
      environment <- environment[sel]
      interventions <- interventions[sel]
    }
   
    # run method
    switch(method,
           "ICP" = {
             result <- runICP(X, environment, interventions, parentsOf, alpha, 
                              variableSelMat, excludeTargetInterventions, 
                              pointConf, setOptions, directed, verbose, ...)
           },
           
           "hiddenICP" = {
             result <- runHiddenICP(X, environment, interventions, parentsOf, 
                                    alpha, variableSelMat, 
                                    excludeTargetInterventions, pointConf, 
                                    setOptions, directed, verbose, ...)
            },
           
            "backShift" = {
              result <- runBackShift(X, environment, parentsOf, 
                                     pointConf, setOptions, verbose, ...)
            },
           
            "regression" = {
              result <- runRegression(X, parentsOf, variableSelMat, pointConf, 
                                      setOptions, directed, verbose, ...)
            },
           
            "gies" = {
              result <- runGIES(X, interventions, parentsOf, variableSelMat, 
                                setOptions, directed, verbose, ...)
            },
           
            "rankGies" = {
              result <- runNonparanormalGIES(X, interventions, parentsOf, 
                                             variableSelMat, setOptions, 
                                             directed, verbose, ...)
            },
            
            "ges" = {
              result <- runGES(X, parentsOf, variableSelMat, setOptions, 
                               directed, verbose, ...)
            },
           
            "rankGes" = {
              result <- runNonparanormalGES(X, parentsOf, variableSelMat, setOptions, 
                              directed, verbose, ...)
            },
           
            "arges" = {
              result <- runARGES(X, parentsOf, variableSelMat, setOptions, 
                              directed, verbose, ...)
            },
           
            "rankArges" = {
              result <- runNonparanormalARGES(X, parentsOf, variableSelMat, setOptions, 
                                directed, verbose, ...)
            },
            
            "pc" = {
              result <- runPC(X, parentsOf, alpha, 
                              variableSelMat, setOptions, 
                              directed, verbose, ...)
            },
           
            "rankPc" = {
              result <- runNonparanormalPC(X, parentsOf, alpha, 
                             variableSelMat, setOptions, 
                             directed, verbose, ...)
            },
           
            "fci" = {
              result <- runFCI(X, suffStat = NULL, parentsOf, alpha, 
                               variableSelMat, setOptions, 
                               directed, verbose, ...)
            },
           
            "rankFci" = {
              result <- runNonparanormalFCI(X, parentsOf, alpha, 
                                          variableSelMat, setOptions, 
                                          directed, verbose, ...)
            },
           
            "rfci" = {
              result <- runRFCI(X, parentsOf, alpha, variableSelMat, setOptions, 
                                directed, verbose, ...)
            },
           
            "fciplus" = {
              result <- runFCIPlus(X, parentsOf, alpha, setOptions, 
                               directed, verbose, ...)
            },
           
           # "directLINGAM" = {
           #   result <- runDirectLINGAM(X, parentsOf, pointConf, 
           #                       setOptions,  
           #                       verbose, ...)
           # },
           
            "LINGAM" = {
              result <- runLINGAM(X, parentsOf, pointConf, setOptions,  
                                  verbose, ...)
            },
           
           "CAM" = {
              result <- runCAM(X, interventions, parentsOf, setOptions, 
                               verbose, ...)
            },
           
            "bivariateCAM" = {
              result <- runBivariateCAM(X, parentsOf, variableSelMat, pointConf,
                                        verbose, ...)
            },
           
            "bivariateANM" = {
              result <- runBivariateANM(X, parentsOf, variableSelMat, pointConf, 
                                        verbose, ...)
            },
           
           "mmhc" = {
              result <- runMMHC(X, parentsOf, alpha, variableSelMat, 
                                setOptions, verbose, 
                                ...)
           },
           
           "RESIT" = {
             result <- runRESIT(X, parentsOf, alpha, setOptions, verbose, 
                               ...)
           },
           
           {
               warning(paste("method ", method," not implemented"))
           }
           )
    
    # bring into correct mode
    result <- changeMode(mode, method, result, length(parentsOf))
   
    # prepare output
    if(returnAsList){
        out <- result$resList 
    }else{
        out <- result$resMat 
      
        if(is.null(out)){
          result <- result$resList
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
          
          resmat <- sparseMatrix(i=rowind,
                                 j=colind,
                                 x=x,
                                 dims=c(ncol(X), length(parentsOf)))
          
          colnames(resmat) <- parentsOf
          
          out <- resmat
          
          if(!sparse){
            out <- as(out, "matrix")
          }
        
        }else if(sparse){
          # cast to sparse matrix
          out <- Matrix(out, sparse = TRUE)
        }
        
        if(is.null(rownames(out)))
          rownames(out) <- if(is.null(colnames(X))) as.character(1:ncol(X)) else colnames(X)
        if(is.null(colnames(out)))
          colnames(out) <- rownames(out)[parentsOf]
        
    }
    
    out
}