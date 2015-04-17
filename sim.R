


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
    if(is.null(layout)) layout <- if(typ=="fg") layout.kamada.kawai(G) else layout.circle(G)
    plot(G, layout=layout, vertex.label=labels,vertex.shape="circle",vertex.label.cex=1.5, vertex.label.color=rgb(0.8,0.1,0.1,0.7),vertex.color="white",vertex.frame.color=rgb(0.8,0.1,0.1,0.5), main=main, edge.color=rgb(0.1,0.1,0.1,0.5), edge.arrow.size=0.7,edge.arrow.width=2,vertex.size=30,vertex.label.dist=0,vertex.label.degree=-pi/2)
    return(layout)
}
plotGraph(A,main="true graph")

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
    ## Ahat <- getParents(X, environment, interventions=interventions,
    ##                        method=method ,alpha=0.1)
    
    ## Option 2): use a stability selection based estimator
    ## with expected number of false positives bounded by EV=2
    Ahat <- getParentsStable(X, environment,EV=2, interventions=interventions,
                             method=method ,alpha=0.1)


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
methods <- c("hiddenICE","ICP","hiddenICP", "lingam", "pc", "rfci","regression","gies","ges")

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





################## Sachs data

data.path <- "../../../feedback/Data/finseldifffifty"
logtransform <- FALSE
G <- 20
A <- diag(6)*0
A[1,c(2,4,5)] <- A[2,c(1,3,5)] <- A[3,2] <- A[4,1] <- A[5,6] <- A[6,5] <- 1

data.path <- "../../../feedback/Data/sachs"
logtransform <- TRUE
G <- 14
A <- diag(11)*0
A[1,2] <- A[2,c(1,6)] <- A[3,c(4,5,9)] <- A[4,c(3,11)] <- A[5,c(3,4,7)] <- A[6,7] <- A[7,6] <- A[8,c(1,2,6,7,10,11)] <- A[9,c(8,10,11,1,2)] <- 1


X <- matrix(nrow=0,ncol=p)
environment <- numeric(0)
for (lc in 1:G){
    tmp <- read.csv( file=paste(data.path,lc,".csv",sep=""), header=TRUE)
    if(logtransform) tmp <- log(tmp)
    if(lc>1) colnames(tmp) <- colnames(X)
    X <- rbind(X, tmp)
    environment <- c(environment, rep(lc,nrow(tmp)))
}

if(length(grep(data.path,pattern="sachs"))>0){
    colnames(X) <- c("Raf","Mek","PLCg","PIP2","PIP3","Erk","Akt","PKA","PKC","p38","JNK")
}
if(length(grep(data.path,pattern="fifty"))>0){
    colnames(X) <- c("S&P 500","10year T-Note","10year Bund", "Oil", "Gold","Silver")
}
####### apply possible  methods (given in vector 'methods')
####### (using all data pooled for pc/lingam/rfci --
#######    --can be changed with option 'onlyObservationalData=TRUE')
methods <- c("hiddenICE","lingam", "pc", "rfci","regression","ges")

## arrange graphical output into a rectangular grid
sq <- ceiling(sqrt(length(methods)+1))
par(mfrow=c(ceiling((length(methods)+1)/sq),sq))
colnames(A) <- colnames(X)
la <- plotGraph(A,main="true graph")
## loop over all methods and compute and print/plot estimate
for (method in methods){
    cat("\n result for method", method,"  ------  \n" )

    ## Option 1): use this estimator as a point estimate if desired:
    ## Ahat <- getParents(X, environment, interventions=interventions,
    ##                        method=method ,alpha=0.1)
    ## Option 2): use a stability selection based estimator
    ## with expected number of false positives bounded by EV=2
    Ahat <- getParentsStable(X, environment,EV=2,  method=method ,alpha=0.1)


    ## print and plot estimate
    print(Ahat)
    plotGraph(Ahat,main=paste("estimate for method",method),layout=la)
}








