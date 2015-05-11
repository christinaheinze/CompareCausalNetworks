\name{getParents}
\alias{getParents}
\title{
Estimate the connectivity matrix of a causal graph
}
\description{
Estimates the connectivity matrix of a directed causal graph, using various possible methods. Supported methods at the moment are PC, Lingam, GES, GIES, RFCI, CAM, and invariant prediction with or without hidden variables (ICP and hiddenICP) as well as invariant covariance estimation (hiddenICE).
}
\usage{
getParents(X, environment = NULL, interventions = NULL, parentsOf = 1:ncol(X), 
          method = c("hiddenICP", "ICP", "hiddenICE", "pc", "lingam",
                      "ges", "gies", "cam", "rfci", "regression", "bivariateANM", "bivariateCAM")[1], 
          alpha = 0.1, variableSelMat = NULL, excludeTargetInterventions = TRUE, 
          onlyObservationalData = FALSE, indexObservationalData = 1, returnAsList = FALSE, confBound = FALSE, 
          setOptions = list(), warnings = TRUE, directed = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A (nxp)-data matrix with n observations of p variables.
}
  \item{environment}{
An optional vector of length n, where the entry for observation i is an index for the environment in which observation i took place (simplest case entries \code{1} for observational data and entries \code{2} for interventional data of unspecified type). Is required for methods "ICP", "hiddenICP", "hiddenICE".
}
  \item{interventions}{
A optional list of length n. The entry for observation i is a numeric vector that specifies the variables on which interventions happened for observation i (a scalar if an intervention happened on just one variable and \code{numeric(0)} if no intervention occured for this observation). Is used for method \code{gies} but will generate the vector \code{environment} if this is set to \code{NULL} (even though it might generate too many different environments for some data so a hand-picked vector \code{environment} is preferable). Is also used for \code{ICP} and \code{hiddenICP} to exclude interventions on the target variable of interest.
}
  \item{parentsOf}{
The variables for which we would like to estimate the parents. Default are all variables.
}
  \item{method}{
 A string that specfies the method to use. The methods \code{pc} (PC-algorithm), \code{lingam} (Lingam), \code{ges} (Greedy equivalence search), \code{gies} (Greedy interventional equivalence search) and \code{rfci} (Really fast causal inference) are imported from the package "pcalg" and are documented there in more detail, including the additional options that can be supplied via \code{setOptions}. The method \code{cam} (Causal additive models) is documented in the package "cam" and the methods \code{ICP} (Invariant causal prediction), \code{hiddenICP} (Invariant causal prediction with hidden variables) and \code{hiddenICE} (Invariant covariance estimation with hidden variables) are from the package "InvariantCausalPrediction". Finally, the methods \code{bivariateANM} and \code{bivariateCAM} are for now implemented internally but will hopefully be part of another package at some point in the near future.
}
  \item{alpha}{
The level at which tests are done. This leads to confidence intervals for \code{ICP} and \code{hiddenICP} and is used internally for \code{pc} and \code{rfci}.
}
  \item{variableSelMat}{
An optional logical matrix of dimension (pxp). An entry \code{TRUE} for entry (i,j) says that variable i should be considered as a potential parent for variable j and vice versa for \code{FALSE}. If the default value of \code{NULL} is used, all variables will be considered, but this can be very slow, especially for methods "pc", "ges", "gies", "rfci" and "cam".
}
  \item{excludeTargetInterventions}{
When looking for parents of variable k in 1,...,p, set to \code{TRUE} if observations where an intervention on variable k occured should be excluded. Default is \code{TRUE}.
}
  \item{onlyObservationalData}{
If set to \code{TRUE}, only observational data is used. It will take the index in \code{environment} specified by \code{indexObservationalData}. If \code{environment} is \code{NULL}, all observations are used. Default is \code{FALSE}.
%If \code{interventions} is specified, it will read out observational data from this list, otherwise it will take the first occuring entry in \code{environment} as a code for observational data. Default is \code{FALSE}.
}
  \item{indexObservationalData}{
Index in \code{environment} that encodes observational data. Default is \code{1}.
}
  \item{returnAsList}{
If set to \code{TRUE}, will return a list, where entry k is a list containing the estimated parents of variable k. The option \code{directed} will be ignored if set to \code{TRUE}. Default is \code{FALSE}.
}
  \item{confBound}{
If \code{TRUE}, numerical estimates will be returned if possible. For methods \code{ICP} and \code{hiddenICP}, these are the values in the individual confidence intervals (at chosen level \code{alpha}) that are closest to 0; for other methods these are point estimates. Some methods do not return numerical point estimates; for these the output will remain binary 0/1 (no-edge/edge). Default is \code{FALSE}.
}
  \item{setOptions}{
A list that can take method-specific options; see the individual documentations of the methods for more options and their possible values.
}
  \item{warnings}{
If set to \code{TRUE}, warnings will be printed. 
}
  \item{directed}{
If \code{TRUE}, an edge will be returned if and only if an edge has been detected to be directed (ie entry will be set to 0 for entry (j,k) if both j->k and k-> j are estimated). Ignored if not the whole graph is estimated or if \code{returnAsList} is \code{TRUE}.
}
}
\value{
If option \code{returnAsList} is \code{FALSE}, a sparse matrix, where a 0 entry in position (j,k) corresponds to an estimate of "no edge" \code{j} -> \code{parentsOf[k]}, while an entry 1 corresponds to an estimated egde. If option \code{confBound} is \code{TRUE}, the 1 entries will be replaced by numerical values that are either point estimates of the causal coefficients or confidence bounds (see above). 

If option \code{returnAsList} is \code{TRUE}, a list will be returned. The k-th entry in the list is the numeric vector with the indices of the estimated parents of node \code{parentsOf[k]}. 
}
\references{
None yet...
}
\author{
Nicolai Meinshausen <meinshausen@stat.math.ethz.ch>
}
\seealso{
\code{\link{getParentsStable}} for stability selection-based estimation of the causal graph. \code{\link{plotGraph}} for plotting the estimated adjacency matrix if the whole graph is estimated. 
}
\examples{
##########################################
######## 1st example:
######## Simulate data with connectivity matrix A with assumptions 
######## 1) hidden variables present
######## 2) precise location of interventions is assumed unknown
######## 3) different environments can be distinguished
##########################################
####### simulate data
set.seed(1)
## sample size n
n <- 1000
## p=3 predictor variables and connectivity matrix A
p  <- 3
A <- diag(p)*0
A[1,2] <- 0.8
A[2,3] <- 0.8
A[3,1] <- -0.4  

## divide data in 10 different environments
G <- 10
environment <- rep(1:G, each=ceiling(n/G))[1:n]
X <- Perturb <-  matrix(0,nrow=n,ncol=p)
## Input of hidden variables into each variable
gamma <- rnorm(p)
Input <- outer(W <- rnorm(n),gamma,FUN="*")
## simulate noise perturbations in each environment
for (i in unique(environment)){
  ind <- which(environment==i)
  multiplier <- rexp(p)*3
  Perturb[ind,] <- sweep(matrix(rnorm(length(ind)*p),ncol=p),
                          2, multiplier,FUN="*")
}
## iterate model to get stable solution
## (necessary only if feedbacks are included)
niter <- 100
for (iter in 1:niter){
  X <- X \%*\% A + Input + Perturb 
}

####### apply all  methods given in vector 'methods'
####### (using all data pooled for pc/lingam/rfci --
#######  -- can be changed with option 'onlyObservationalData=TRUE')
methods <- c("hiddenICE", "lingam", "pc", "rfci","regression")

## arrange graphical output into a rectangular grid
sq <- ceiling(sqrt(length(methods)+1))
par(mfrow=c(ceiling((length(methods)+1)/sq),sq))

## plot and print true graph
cat("\n true graph is  ------  \n" )
print(A)
plotGraph(A,main="true graph")

## loop over all methods and compute and print/plot estimate
for (method in methods){
  cat("\n result for method", method,"  ------  \n" )
  
  ## Option 1): use this estimator as a point estimate
#   Ahat <- getParents(X, environment, method=method, alpha=0.1)
  
  ## Option 2): use a stability selection based estimator
  ## with expected number of false positives bounded by EV=2
  Ahat <- getParentsStable(X, environment, EV=2, method=method ,alpha=0.1)
  
  ## print and plot estimate
  print(Ahat)
  plotGraph(Ahat,main=paste("estimate for method",method))
}





##########################################
######## 2nd example:
######## Simulate data with connectivity matrix A with assumptions
######## 1) No hidden variables
######## 2) Precise location of interventions is known
##########################################
####### simulate data
set.seed(1)
## sample size n
n <- 2000
## p=5 predictor variables
p  <- 5
A <- diag(p)*0
A[1,2] <- 0.8
A[2,3] <- -0.8
A[3,4] <- 0.8
A[3,5] <- 0.8
A[4,5] <- 0.3
## can add/remove feedback by using/not using
A[5,2] <- 0.8 

## choose explicity intervention targets
interventions <- list()
for (i in 1:n) interventions[[i]] <- sample(1:p,2)
environment <- match(interventions, unique(interventions))
X <- Perturb <-  matrix(0,nrow=n,ncol=p)
## Independent noise at each variable 
Input <- matrix( rnorm(n*p),nrow=n)
## change level of noise for each intervention
for (i in 1:n){
    Perturb[i, interventions[[i]]] <- rnorm(length(interventions[[i]]))*5
}
## iterate model to get stable solution 
## (necessary only if feedbacks are included)
niter <- 100
for (iter in 1:niter){
    X <- X \%*\% A + Input + Perturb 
}


####### apply possible  methods given in vector 'methods'
####### (using all data pooled for pc/lingam/rfci --
#######    --can be changed with option 'onlyObservationalData=TRUE')
methods <- c("hiddenICE","ICP","hiddenICP", "lingam", "pc", "rfci",
              "regression","gies","ges")

## arrange graphical output into a rectangular grid
sq <- ceiling(sqrt(length(methods)+1))
par(mfrow=c(ceiling((length(methods)+1)/sq),sq))

## plot and print true graph
cat("\n true graph is  ------  \n" )
print(A)
plotGraph(A,main="true graph")

## loop over all methods and compute and print/plot estimate
for (method in methods){
    cat("\n result for method", method,"  ------  \n" )

    ## Option 1): use this estimator as a point estimate if desired:
#     Ahat <- getParents(X, environment, interventions=interventions,
#                       method=method ,alpha=0.1)

    ## Option 2): use a stability selection based estimator
    ## with expected number of false positives bounded by EV=2
    Ahat <- getParentsStable(X, environment,EV=2, interventions=interventions, 
                            method=method ,alpha=0.1)


    ## print and plot estimate
    print(Ahat)
    plotGraph(Ahat,main=paste("estimate for method",method))
}



}
\keyword{ Causality }
\keyword{ Graph estimation }
