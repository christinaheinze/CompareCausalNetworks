runNonparanormalARGES <- function(X, parentsOf, variableSelMat, setOptions, directed, verbose, 
                     result, ...){
  
  package <- setOptions$package
  if(is.null(package)) package <- "huge"
  
  if(package == "huge"){
    method <- setOptions$method
    if(is.null(method)) method <- "mb"
    criterion <- setOptions$criterion
    if(is.null(criterion)) criterion <- "ric"
  }else if(package == "flare"){
    method <- setOptions$method
    if(is.null(method)) method <- "tiger"
    criterion <- setOptions$criterion
    if(is.null(criterion)) criterion <- "cv"
  }else{
    stop(paste("Package", package, "not supported for CIG estimation. Valid
               options are 'huge' or 'flare'."))
  }
  
  # estimate CIG
  if(is.null(variableSelMat)){
    if(package == "huge"){
      hugeObj <- huge::huge(X, method = method, verbose = FALSE)
      hugeSel <- huge::huge.select(hugeObj, criterion = criterion, verbose = FALSE)
      variableSelMat <- hugeSel$refit
    }else{
      flareObj <- flare::sugm(X, method = method, verbose = FALSE)  
      flareSel <- flare::select(flareObj,  criterion = criterion, verbose = FALSE)  
      variableSelMat <- flareSel$refit
    }
    
    variableSelMat <- as.matrix(variableSelMat)
    variableSelMat[variableSelMat == 1] <- TRUE
    variableSelMat[variableSelMat == 0] <- FALSE
  }
  
  # additional options for ARGES
  if(is.null(setOptions$adaptive)){
    setOptions$adaptive <- "vstructures" #ARGES-CIG
  }else{
    if(setOptions$adaptive == "none")
      setOptions$adaptive <- "vstructures" #ARGES-CIG
  }
  
  given.cov.mat <- cov(X)
  p <- ncol(X)
  n <- nrow(X)
  e <- eigen(given.cov.mat)
  sqrt.given.cov.mat <- e$vectors%*%sqrt(diag(e$values))
  
  ## generate from N(0,I)
  dat <- rmvnorm(n,mean=rep(0,p),sigma=diag(p))
  
  ## transform data so that cov(dat) = given.cov.mat
  samp.cov.mat <- 2*sin(cor(dat,method="spearman")*pi/6)
  e <- eigen(samp.cov.mat, symmetric = T)
  e$values[which(e$values<0)] <- 0
  samp.cov.mat <- e$vectors%*%diag(e$values)%*%t(e$vectors)
  e <- eigen(samp.cov.mat)
  sqrt.samp.cov.mat <- e$vectors%*%sqrt(diag(e$values))
  X <- t(sqrt.given.cov.mat%*%solve(sqrt.samp.cov.mat,t(dat)))
  
  runGES(X, parentsOf, variableSelMat, setOptions, directed, verbose, 
         result, ...)
}