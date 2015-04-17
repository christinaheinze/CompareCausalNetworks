


library(InvariantCausalPrediction)
library(pcalg)
require(clue)
require(LINGAM)
require(CAM)
require(linprog)
require(glmnet)
require(igraph)
require(jointDiag)
source("getParents.R")
source("getParentsStable.R")
source("computeDelta.R")
source("edgeRetention.R")
source("edgeSelection.R")
source("bivariateCAM.R")
source("bivariateANM.R")
source("train_gam.R")
source("indtestHsic.R")

lf <- list.files("../InvariantCausalPrediction/R")
lf <- lf[grep(lf,pattern="*.R$")]
for (k in lf) source(paste("../InvariantCausalPrediction/R/",k,sep=""))
##########################################
######## 1st example:
######## Simulate data with connectivity matrix A with assumptions 
######## 1) hidden variables present
######## 2) precise location of interventions is assumed unknown
######## 3) but different environments can be distinguished
##########################################


plotGraph <- function(A, main="",labels=NULL,layout=NULL,typ="fg"){
    if(is.null(labels)) labels <- if( !is.null(cc <- colnames(A)))
        cc else as.character(1:ncol(A))
    G <- graph.adjacency(A,mode="directed",weighted="a")
    plot(G, layout=layout.circle(G), vertex.label=labels,vertex.shape="circle",vertex.label.cex=1.5, vertex.label.color=rgb(0.8,0.1,0.1,0.7),vertex.color="white",vertex.frame.color=rgb(0.8,0.1,0.1,0.5), main=main, edge.color=rgb(0.1,0.1,0.1,0.5), edge.arrow.size=0.7,edge.arrow.width=2,vertex.size=30,vertex.label.dist=0,vertex.label.degree=-pi/2)
}


set.seed(1)
## sample size n
n <- 1000
## p=3 predictor variables and connectivity matrix A
p  <- 3
A <- diag(p)*0
A[1,2] <- 0.8
A[2,3] <- 0.8
A[3,1] <- -0.4

## divide data into 10 different environments
G <- 10
environment <- rep(1:G, each=ceiling(n/G))[1:n]
X <- Perturb <-  matrix(0,nrow=n,ncol=p)
## Input of hidden variables into each variable
alpha <- rnorm(p)
Input <- outer(W <- rnorm(n),alpha,FUN="*")
## simulate noise perturbations in each environment
for (i in unique(environment)){
    ind <- which(environment==i)
    multiplier <- rexp(p)*3
    Perturb[ind,] <- sweep( matrix(rnorm(length(ind)*p),ncol=p),
                           2, multiplier,FUN="*")
}
## iterate model to get stable solution
niter <- 100
for (iter in 1:niter){
    X <- X %*% A + Input + Perturb 
}

####### apply all possible  methods (given in vector 'methods')
####### (using all data pooled for pc/lingam/rfci --
#######  -- can be changed with option 'onlyObservationalData=TRUE')

methods <- c("hiddenICE", "lingam", "pc", "rfci","regression","bivariateANM","bivariateCAM")

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
    ## Ahat <- getParents(X, environment, interventions=interventions,
    ##                        method=method ,alpha=0.1)
    
    ## Option 2): use a stability selection based estimator
    ## with expected number of false positives bounded by EV=2
    Ahat <- getParentsStable(X, environment,EV=2, method=method ,alpha=0.1)


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
## iterate model to get stable solution (necessary only if feedbacks are included)
niter <- 100
for (iter in 1:niter){
    X <- X %*% A + Input + Perturb 
}


####### apply possible  methods (given in vector 'methods')
####### (using all data pooled for pc/lingam/rfci --
#######    --can be changed with option 'onlyObservationalData=TRUE')
methods <- c("hiddenICE","ICP","hiddenICP", "lingam", "pc", "rfci","regression","gies","ges","bivariateCAM","bivariateANM")

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
    ## Ahat <- getParents(X, environment, interventions=interventions,
    ##                        method=method ,alpha=0.1)
    ## Option 2): use a stability selection based estimator
    ## with expected number of false positives bounded by EV=2
    Ahat <- getParentsStable(X, environment,EV=2, interventions=interventions, method=method ,alpha=0.1)


    ## print and plot estimate
    print(Ahat)
    plotGraph(Ahat,main=paste("estimate for method",method))
}




