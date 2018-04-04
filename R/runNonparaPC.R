runNonparanormalPC <- function(X, parentsOf, alpha, variableSelMat, setOptions, directed, verbose, 
                               ...){
  
  suffStat.data <- list(C=2*sin(cor(X,method="spearman")*pi/6), n=nrow(X))
  setOptions$indepTest <- pcalg::gaussCItest
  setOptions$u2pd <- "relaxed"
  setOptions$suffStat <- suffStat.data
  
  runPC(X, parentsOf, alpha, variableSelMat, setOptions, 
        directed, verbose, ...)
}