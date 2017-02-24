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
                                  sparse, numberInt, strengthInt,
                                  cyclic, strengthCycle,
                                  seed =1){
###### set seed
set.seed(seed)
  
  # generate A
  cont <- TRUE
  while(cont){
    ## skeleton
    AS <- 0*diag(p)
    for (k in 1:p){
      for (k2 in 1:p){
        if(k2>k) AS[k,k2] <- rbinom(1,1,sparse)
        # if(k2<k) AS[k2,k] <- rbinom(1,1,sparseBack)
      }
    }
    
    Aorig <- AS * matrix( runif(p^2,-1,1),nrow=p)
    # eig <- max(abs(eigen(Aorig)$values))
    # if(eig>0.5) Aorig <- Aorig/(0.5+eig)
    if(sum(Aorig!=0)>0) cont <- FALSE
  }
  A <- Aorig
  # print(A)

  
  ####### initialize
  Perturb <- matrix(0,nrow=n,ncol=p)

  ### simulate environments
  q <- floor(p/numberInt)
  whereInt <- list()
  allVar <- 1:p
  for (nic in 1:numberInt){
    if(nic==numberInt){
      whereInt[[nic]] <- allVar
    }else{
      whereInt[[nic]] <- sample(allVar,q)  
      allVar <- allVar[ !(allVar %in% whereInt[[nic]])]
    }
  }
  
  environment <- sample(0:numberInt, 
                        n, prob=c(1/2, rep(1/(2*numberInt), numberInt)),
                        replace=TRUE)
  interventions <- list()
  for (i in 1:n){
    interventions[[i]] <- numeric(0)
    if(environment[i]>0) interventions[[i]] <- whereInt[[environment[i]]]
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
  A <- setSNR(noise = noise[environment == 0,], A)
    
  # print(A)
  if(cyclic){
    # add cycle to A
    # sample uniformly at random one connection backwards
    backConnect <- sample(1:p, 2, replace = FALSE)
    A[max(backConnect), min(backConnect)] <- 1
    # find cyclic path from max(backConnect) back to max(backConnect)
    # and get largest coefficient path can have (if it exists)
    cpmat <- getLargestWeightForCycle(A)$cpMat
    # sample sign of weight of backward connection
    signEntry <- sample(c(-1,1),1)
    # adjust edge weight by strengthCycle parameter
    if(cpmat[max(backConnect), max(backConnect)] != 0){
      A[max(backConnect), min(backConnect)] <- signEntry*strengthCycle/
        cpmat[max(backConnect), max(backConnect)]
    }else{
      A[max(backConnect), min(backConnect)] <- runif(1,-1,1)
    }
   # print(A)
    
  }
  
  # apply t-distribution for noise
  noise <- apply(noise, 2, function(x) qt(pnorm(x), df = df))
  inv <- solve(diag(p) - A)
  X <- (noise + Perturb)%*%inv

  config.options <- list(trueA = A, 
                         n = n, 
                         p = p,
                         df = df, 
                         rhoNoise = rhoNoise, 
                         snrPar = snrPar,
                         sparse = sparse,
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