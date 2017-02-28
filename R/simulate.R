#' Simulate data of a causal cyclic model under shift interventions.
#' 
#' @description Simulate data of a causal cyclic model under shift interventions.
#' 
#' @param n Number of observations.
#' @param p Number of variables.
#' @param df
#' @param rhoNoise
#' @param snrPar
#' @param sparse
#' @param numberInt 
#' @param strengthInt Regulates the strength of the interventions.
#' @param cyclic
#' @param strengthCycle 
#' @param seed Random seed.
#' @return A list with the following elements: 
#' \itemize{
#'   \item \code{X} (nxp)-dimensional data matrix
#'   \item \code{environment} Indicator of the experiment or the intervention type an 
#'   observation belongs to. A numeric vector of length n. 
#'   \item \code{interventions} Location of interventions
#'   \item \code{whereInt}
#'   \item \code{noise}
#'   \item \code{configs} A list with the following elements: 
#'   \itemize{
#'   \item ...
#'   }
#' }
simulateInterventions <- function(n, p, df, rhoNoise, snrPar,
                                  sparse, doInterv,
                                  numberInt, strengthInt,
                                  cyclic, strengthCycle,
                                  seed =1){
  ###### set seed
  set.seed(seed)
  
  if(snrPar == 0 | snrPar > 1){
    stop("snrPar needs to be larger than 0 and smaller or equal to 1.")
  }
  
  if(strengthCycle >= 1){
    stop("strengthCycle needs to be smaller than 1.")
  }
  
  # generate A
  ## skeleton
  AS <- 0*diag(p)
  for (k in 1:p){
    for (k2 in 1:p){
      if(k2>k) AS[k,k2] <- rbinom(1,1,sparse)
    }
  }
    
  A <- AS * matrix( runif(p^2,-1,1),nrow=p)
  
  
  ####### initialize
  Perturb <- matrix(0,nrow=n,ncol=p)

  ### simulate environments
  environment_var <- sample(1:numberInt, p, replace=TRUE)
  
  whereInt <- list()
  for(nic in 1:numberInt){
    whereInt[[nic]] <- which(environment_var == nic)  
  }

  environment <- sample(1:numberInt, n, replace=TRUE)
  interventions <- list()
  for (i in 1:n){
    interventions[[i]] <- whereInt[[environment[i]]]
  }
  
  ###### simulate shift interventions 
  for(intervEnv in 1:numberInt){
    rowsSelect <- environment == intervEnv
    colsSelect <- whereInt[[intervEnv]]
    Perturb[rowsSelect,  whereInt[[intervEnv]]] <-  
      strengthInt*rt( length(colsSelect) * sum(rowsSelect) ,df=df)
  }
  
  ###### simulate noise
  noise <- matrix(rnorm(n*p),nrow=n)
  Sigma <- matrix(rhoNoise,nrow=p,ncol=p)
  diag(Sigma) <- 1
  CC <- chol(Sigma)
  noise <- noise %*% CC
  noise[,which(colSums(A) != 0)] <- sqrt(snrPar)*noise[,which(colSums(A) != 0)] 
  
  # change A to achieve given SNR in observational case
  A <- setSNR(noise = noise, A)
    
  if(cyclic){
    cont <- TRUE
    Atmp <- A
    while(cont){
      # add cycle to A
      # sample uniformly at random one connection backwards
      backConnect <- sample(1:p, 2, replace = FALSE)
      Atmp[max(backConnect), min(backConnect)] <- 1
      # find cyclic path from max(backConnect) back to max(backConnect)
      # and get largest coefficient path can have (if it exists)
      cpmat <- getLargestWeightForCycle(Atmp)$cpMat
      if(cpmat[max(backConnect), max(backConnect)] != 0){
        cont <- FALSE
      }else{
        Atmp <- A
      }
    }
    
    # sample sign of weight of backward connection
    signEntry <- sample(c(-1,1),1)
    # adjust edge weight by strengthCycle parameter
    A[max(backConnect), min(backConnect)] <- signEntry*strengthCycle/
        cpmat[max(backConnect), max(backConnect)]
  }
  
  # apply t-distribution for noise
  noise <- apply(noise, 2, function(x) qt(pnorm(x), df = df))
  
  if(doInterv){
    X <- matrix(0, nrow = n, ncol = p)
    # for do interventions cycle through environments and adjust A 
    # as well as noise + Perturb accordingly
    for(env in 1:numberInt){
      # data points from setting env
      inds <- which(environment == env) 
      # modify A: intervention on var j means to cut all incoming connections
      Amod <- A
      Amod[,whereInt[[env]]] <- 0
      # get inverse
      inv <- solve(diag(p) - Amod) 
      # replace columns of noise where intervention was applied
      noiseMat <- noise[inds,,drop = FALSE]
      noiseMat[,whereInt[[env]]] <- Perturb[inds,whereInt[[env]], drop = FALSE]
      X[inds,] <- noiseMat%*%inv
    }     
  }else{
    inv <- solve(diag(p) - A)
    X <- (noise + Perturb)%*%inv
  }
  
  config.options <- list(trueA = A, 
                         n = n, 
                         p = p,
                         df = df, 
                         rhoNoise = rhoNoise, 
                         snrPar = snrPar,
                         sparse = sparse,
                         doInterv = doInterv,
                         numberInt = numberInt,
                         strengthInt = strengthInt,
                         cyclic = cyclic, 
                         strengthCycle = strengthCycle,
                         seed = seed)
                         

  list(X = X, 
       environment = environment,
       interventions = interventions,
       whereInt = whereInt,
       noise = noise,
       configs = config.options)
}