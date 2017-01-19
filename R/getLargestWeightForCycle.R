getLargestWeightForCycle <- function(AdjMat, verbose = FALSE){
  cp.mat <- abs(AdjMat)
  p.relevant <- ncol(AdjMat)
  nodes.in.path.from.i.to.i <- vector("list", length = p.relevant)
  
  for(k in 1:p.relevant){
    for(i in 1:p.relevant){
      for(j in 1:p.relevant){
        path.over.k <- cp.mat[i,k] * cp.mat[k,j]
        path.not.over.k <- cp.mat[i,j]
        
        cp.mat[i,j] <- max(path.not.over.k, path.over.k)
        if(i == j){
          
          if(path.over.k > path.not.over.k){
            nodes.in.path.from.i.to.i[[i]] <- c(nodes.in.path.from.i.to.i[[i]], k)
            if(verbose){
              cat('Going from', i, 'to', j, 'over', k, '.\n')
            }
          }
          
          # if now the cycle product is >= 1, found inadmissible model
          if(cp.mat[i,j] >= 1){
            if(verbose){
              cat("Cycle found with length", cp.mat[i,j], "at node", i, 
                  ". Model is not stable. \n")
              cat('Involved inner nodes are', 
                  nodes.in.path.from.i.to.i[[i]], '.\n')
            }
            nodes.involved.in.cycle <- c(i, nodes.in.path.from.i.to.i[[i]])
            return(list(success = FALSE, cycleNodes = nodes.involved.in.cycle, cpMat = cp.mat))
          }
        }
      }
    }
  }
  list(success = TRUE, cycleNodes = NULL, cpMat = cp.mat)
}