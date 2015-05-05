
edgeRetention <- function(matrixAList, threshold, p){
  summedAdjacency <- matrix(0, p, p)

  Ahata <- array( unlist(matrixAList), dim=c( nrow(matrixAList[[1]]), ncol(matrixAList[[1]]), length(matrixAList)))
  Ahatpos <- apply(Ahata, c(1,2), function(x) mean(x!=0))
  summedAdjacency <- Ahatpos

  summedAdjacency[summedAdjacency < threshold] <- 0
  
  round(100*summedAdjacency)
}
