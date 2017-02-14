convertForRanking <- function(amat, queries, method){
  # queriesList <- vector("list", length = length(queries))
  # names(queriesList) <- queries
  
  # get graph mode to know how to interpret encoding
  graphMode <- getGraphMode(method)
  
  if(graphMode == "PAG"){
    queriesList <- lapply(queries, function(q) answerQueries(amat, NULL, q))
    
  }else if(graphMode == "CPDAG"){
    # to answer ancestral queries -- convert to ancestral matrix
    ancMat <- cpdag2pag(amat)
    queriesList <- lapply(queries, function(q) answerQueries(ancMat, amat, q))
    
  }else if(graphMode == "DAG"){
    # to answer ancestral queries -- convert to ancestral matrix
    ancMat <- dag2ancestral(amat)
    queriesList <- lapply(queries, function(q) answerQueries(ancMat, amat, q))
    
  }else if(graphMode == "CG"){
    # to answer ancestral queries -- convert to ancestral matrix
    ancMat <- cg2ancestral(amat)
    queriesList <- lapply(queries, function(q) answerQueries(ancMat, amat, q))
  }
  
  # check query
  names(queriesList) <- queries
  
  # return: four matrices with increment to hist
  queriesList
}

answerQueries <- function(ancestralAmat, parentalAmat, query){
  p <- ncol(ancestralAmat)
  
  if(graphMode == "PAG"){
    # can answer ancestral queries
    switch(query,
           "isParent" = { 
             matrix(0, p, p)
           },
           "isMaybeParent" = { 
             matrix(0, p, p)
           }, 
           "isNoParent" = { 
             matrix(0, p, p)
           },
           "isAncestor" = { 
             resMat <- pagIsAncestor(ancestralAmat)
           },
           "isMaybeAncestor" = { 
             resMat <- pagIsMaybeAncestor(ancestralAmat)
            }, 
           "isNoAncestor" = { 
             resMat <- pagIsNoAncestor(ancestralAmat)
            })
  }else if(graphMode == "CPDAG"){
    # to answer ancestral queries -- convert to ancestral matrix
    switch(query,
           "isParent" = { 
             resMat <- parentalAmat * (t(parentalAmat)==0)
             resMat[resMat != 0] <- 1
           },
           "isMaybeParent" = { 
             resMat <- parentalAmat * (t(parentalAmat)!=0)
             resMat[resMat != 0] <- 1
           }, 
           "isNoParent" = { 
             noParentMat <- matrix(0, ncol = p, nrow = p)
             defParentalMat <- parentalAmat * (t(parentalAmat)==0)
             maybeParentalMat <- parentalAmat * (t(parentalAmat)!=0)
             noParentMat[defParentalMat == 0] <- 1
             noParentMat[maybeParentalMat == 1] <- 0
             resMat <- noParentMat
           },
           "isAncestor" = { 
             resMat <- pagIsAncestor(ancestralAmat)
           },
           "isMaybeAncestor" = { 
             resMat <- pagIsMaybeAncestor(ancestralAmat)
           }, 
           "isNoAncestor" = { 
             resMat <- pagIsNoAncestor(ancestralAmat)
           })
    
  }else if(graphMode == "DAG" | graphMode == "CG"){
    # to answer ancestral queries -- convert to ancestral matrix
    switch(query,
           "isParent" = { 
             resMat <- parentalAmat 
             resMat[resMat != 0] <- 1
           },
           "isMaybeParent" = { 
             matrix(0, p, p)
           }, 
           "isNoParent" = { 
             noParentMat <- matrix(0, ncol = p, nrow = p)
             noParentMat[parentalAmat == 0] <- 1
             resMat <- noParentMat
           },
           "isAncestor" = { 
             resMat <- ancestralAmat
           },
           "isMaybeAncestor" = { 
             matrix(0, p, p)
           }, 
           "isNoAncestor" = { 
             noParentMat <- matrix(0, ncol = p, nrow = p)
             noParentMat[ancestralAmat == 0] <- 1
             resMat <- noParentMat
           })
  }else if(graphMode == "CG"){
    # to answer ancestral queries -- convert to ancestral matrix
    switch(query,
           "isParent" = { 
             resMat <- parentalAmat 
             resMat[resMat != 0] <- 1
           },
           "isMaybeParent" = { 
             matrix(0, p, p)
           }, 
           "isNoParent" = { 
             noParentMat <- matrix(0, ncol = p, nrow = p)
             noParentMat[parentalAmat == 0] <- 1
             resMat <- noParentMat
           },
           "isAncestor" = { 
             resMat <- ancestralAmat
           },
           "isMaybeAncestor" = { 
             matrix(0, p, p)
           }, 
           "isNoAncestor" = { 
             noParentMat <- matrix(0, ncol = p, nrow = p)
             noParentMat[ancestralAmat == 0] <- 1
             resMat <- noParentMat
           })
  }
  resMat
}

# dag2ancestral, cyclic2ancestral 
# convert dag to ancestral graph: is there a path along the directed edges in the DAG
# from node i to node j, if so set mat[i,j] = 1 (i.e. i is an ancestor of j)

dag2ancestral <- function(amat){
  g <- cg2ancestral(amat)
  diag(g) <- 0
  g
}

cg2ancestral <- function(amat){
  p <- ncol(amat)
  cumsum <- matrix(0, p,p)
  paths <- lapply(1:p, function(j) amat%^%j)
  for(i in 1:p) cumsum <- cumsum + paths[[i]]
  cumsum[cumsum != 0] <- 1
  cumsum
}

cpdag2pag <- function(amat){
  directedOnly <- amat * (t(amat)==0)
  directedPaths <- dag2ancestral(directedOnly)
  directedAndUndirectedPaths <- dag2ancestral(amat)
  # only directed paths needed --> convert to pag 2-3 edge
  sumMat <- sumMatTmp <- directedPaths + directedAndUndirectedPaths
  sumMatTmp[sumMatTmp != 2] <- 0
  sumMatTmp[t(sumMatTmp) == 2] <- 3
  # combination needed --> convert to maybe edge
  sumMat[sumMat != 1] <- 0
  
  sumMat[sumMat == 1 & t(sumMat) == 0] <- 2
  sumMat[t(sumMat == 2)] <- 1
  
  sumMat + sumMatTmp
}


pagIsAncestor <- function(ancestralAmat){
  resMat <- ancestralAmat
  # find edges of type
  #  amat[a,b] = 2  and  amat[b,a] = 3  implies   a --> b.
  resMatTmp <- resMat + t(resMat)
  resMatTmp[resMatTmp != 5] <- 0
  # remove entry so that resulting resMat[i,j] = 1 
  # means that i is an ancestor of j
  resMatTmp[resMat != 2] <- 0
  resMatTmp[resMatTmp != 0] <- 1
  resMatTmp
}

pagIsMaybeAncestor <- function(ancestralAmat){
  resMat <- ancestralAmat
  # find edges of type
  # amat[a,b] = 1  and  amat[b,a] = 1   implies   a o-o b.
  # amat[a,b] = 1  and  amat[b,a] = 2   implies   a <-o b.
  # amat[a,b] = 1  and  amat[b,a] = 3   implies   a --o b.
  resMatTmp <- resMat + t(resMat)
  resMatTmp[resMat != 1 & t(resMat) != 1] <- 0
  
  # Case 1: a o-o b
  resMatTmpCase1 <- resMatTmp
  resMatTmpCase1[resMatTmpCase1 != 2] <- 0
  
  # Case 2: a <-o b
  resMatTmpCase2 <- resMatTmp
  resMatTmpCase2[resMatTmpCase2 != 3] <- 0
  # remove entry so that resulting resMat[i,j] = 1 
  # means that i is an ancestor of j
  resMatTmpCase2[resMat != 2] <- 0
  
  # Case 3: a --o b
  resMatTmpCase3 <- resMatTmp
  resMatTmpCase3[resMatTmpCase3 != 4] <- 0
  resMatTmpCase3[resMat != 1] <- 0
  
  resMatTmpAll <- resMatTmpCase1 + resMatTmpCase2 + resMatTmpCase3
  resMatTmpAll[resMatTmpAll != 0] <- 1
  resMatTmpAll
}

pagIsNoAncestor <- function(ancestralAmat){
  defAncestralMat <- pagIsAncestor(ancestralAmat)
  maybeAncestralMat <- pagIsMaybeAncestor(ancestralAmat)
  noAncestorMat <- 0*ancestralAmat
  noAncestorMat[defAncestralMat == 0] <- 1
  noAncestorMat[maybeAncestralMat == 1] <- 0
  noAncestorMat
}
