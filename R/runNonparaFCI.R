runNonparanormalFCI <- function(X, parentsOf, alpha, variableSelMat, setOptions, directed, verbose, 
                               ...){
  
  suffStat.data <- list(C=2*sin(cor(X,method="spearman")*pi/6), n=nrow(X))
  setOptions$indepTest <- pcalg::gaussCItest

  runFCI(X, suffStat.data, parentsOf, alpha, variableSelMat, setOptions, 
        directed, verbose, ...)
}