convertForRanking <- function(amat, queries, method, assumeNoSelectionVars = TRUE){
  # queriesList <- vector("list", length = length(queries))
  # names(queriesList) <- queries
  
  # get graph mode to know how to interpret encoding
  graphMode <- getGraphMode(method)
  
  if(graphMode == "PAG"){
    queriesList <- lapply(queries, function(q) 
      answerQueries(amat, NULL, q, graphMode, assumeNoSelectionVars))
    
  }else if(graphMode == "CPDAG"){
    # to answer ancestral queries -- convert to ancestral matrix
    ancMat <- cpdag2pag(amat)
    queriesList <- lapply(queries, function(q) 
      answerQueries(ancMat, amat, q, graphMode, assumeNoSelectionVars))
    
  }else if(graphMode == "DAG"){
    # to answer ancestral queries -- convert to ancestral matrix
    ancMat <- dag2ancestral(amat)
    queriesList <- lapply(queries, function(q) 
      answerQueries(ancMat, amat, q, graphMode, assumeNoSelectionVars))
    
  }else if(graphMode == "CG"){
    # to answer ancestral queries -- convert to ancestral matrix
    ancMat <- cg2ancestral(amat)
    queriesList <- lapply(queries, function(q) 
      answerQueries(ancMat, amat, q, graphMode, assumeNoSelectionVars))
  }
  
  # check query
  names(queriesList) <- queries
  
  # return: four matrices with increment to hist
  queriesList
}

answerQueries <- function(ancestralAmat, parentalAmat, query, graphMode, assumeNoSelectionVars = TRUE){
  if(!is.null(ancestralAmat)) p <- ncol(ancestralAmat) else  p <- ncol(parentalAmat) 
  
  if(graphMode == "PAG"){
    # can answer ancestral queries
    switch(query,
           "isParent" = { 
             resMat <- 0*ancestralAmat
           },
           "isMaybeParent" = { 
             resMat1 <- pagIsAncestorEdge(ancestralAmat, assumeNoSelectionVars)
             resMat2 <- pagIsMaybeAncestorEdge(ancestralAmat, assumeNoSelectionVars)
             resMat <- resMat1 + resMat2
             resMat[resMat != 0] <- 1
           }, 
           "isNoParent" = { 
             maybeParentMat <- pagIsAncestorEdge(ancestralAmat, assumeNoSelectionVars) +
               pagIsMaybeAncestorEdge(ancestralAmat, assumeNoSelectionVars)
             maybeParentMat[maybeParentMat != 0] <- 1
             resMat <- 0*ancestralAmat
             # setting all entries to 1 that are neither parents nor possible parents
             resMat[maybeParentMat == 0] <- 1
           },
           "isAncestor" = { 
             resMat <- pagIsAncestor(ancestralAmat, assumeNoSelectionVars)
           },
           "isMaybeAncestor" = { 
             resMat <- pagIsMaybeAncestor(ancestralAmat, assumeNoSelectionVars)
            }, 
           "isNoAncestor" = { 
             resMat <- pagIsNoAncestor(ancestralAmat, assumeNoSelectionVars)
            })
  }else if(graphMode == "CPDAG"){
    # to answer ancestral queries -- convert to ancestral matrix
    switch(query,
           "isParent" = { 
             resMat <- parentalAmat * (t(parentalAmat)==0)
             resMat[resMat != 0] <- 1
           },
           "isMaybeParent" = { 
             resMat <- parentalAmat 
             resMat[resMat != 0] <- 1
           }, 
           "isNoParent" = { 
             resMat <- 0*parentalAmat
             # setting all entries to 1 that are neither parents nor possible parents
             resMat[parentalAmat == 0] <- 1 
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
             resMat <- parentalAmat
             resMat[resMat != 0] <- 1
           }, 
           "isNoParent" = { 
             resMat <- 0*parentalAmat
             # setting all entries to 1 that are not parents
             resMat[parentalAmat == 0] <- 1
           },
           "isAncestor" = { 
             resMat <- ancestralAmat
           },
           "isMaybeAncestor" = { 
             resMat <- ancestralAmat
           }, 
           "isNoAncestor" = { 
             resMat <- 0*ancestralAmat
             resMat[ancestralAmat == 0] <- 1
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
  # combination of directed and undirected edges needed --> convert to maybe edge 
  sumMat[sumMat != 1] <- 0
  
  sumMat[sumMat == 1 & t(sumMat) == 0] <- 2
  sumMat[t(sumMat == 2)] <- 1
  
  sumMat + sumMatTmp
}

pagIsAncestorEdge <- function(ancestralAmat, assumeNoSelectionVars = TRUE){
  resMat <- ancestralAmat
  
  if(assumeNoSelectionVars){
    # find edges of type --*
    resMat[resMat != 3] <- 0
    resMat <- t(resMat)
    resMat[resMat != 0] <- 1
  }else{
    # find edges of type -->
    resMatTmpCase1 <- resMat + t(resMat)
    resMatTmpCase1[resMatTmpCase1 != 5] <- 0
    # remove entry so that resulting resMat[i,j] = 1
    # means that i is an ancestor of j
    resMatTmpCase1[resMat != 2] <- 0
    resMatTmpCase1[resMatTmpCase1 != 0] <- 1
    resMat <- resMatTmpCase1
  }
  
  resMat
}

pagIsMaybeAncestorEdge <- function(ancestralAmat, assumeNoSelectionVars = TRUE){
  resMat <- ancestralAmat
  
  if(assumeNoSelectionVars){
    # find edges of type o-*
    resMat[resMat != 1] <- 0
    resMat <- t(resMat)
    resMat[resMat != 0] <- 1
    resMat
  }else{
    # find edges of type o-o, o-> and --o
    
    # Case: a o-o b
    resMatTmpCase1 <- resMat + t(resMat)
    # get only edges of type a o-o b
    resMatTmpCase1[resMatTmpCase1 != 2] <- 0

    # Case: a o-> b
    resMatTmpCase2 <- resMat + t(resMat)
    # get only edges of type a o-> b
    resMatTmpCase2[resMatTmpCase2 != 3] <- 0
    # remove entry so that resulting resMat[i,j] = 1
    # means that i might be an ancestor of j
    # -> only keep entries of type 2 so that
    # amat[b,a] = 2 remain --> b might be an ancestor of a (a has edge mark arrowhead,
    # b has edge mark circle)
    resMatTmpCase2[resMat != 2] <- 0
    
    # Case: a --o b
    resMatTmpCase3 <- resMat + t(resMat)
    # retain edges containing circles
    resMatTmpCase3[resMat != 1 & t(resMat) != 1] <- 0
    # get only edges of type  a --o b
    resMatTmpCase3[resMatTmpCase3 != 4] <- 0
    # remove entry so that resulting resMat[i,j] = 1
    # means that i might be an ancestor of j
    # -> only keep entries of type 1 so that
    # amat[a,b] = 1 remain --> a might be an ancestor of b (b has edge mark circle,
    # a has edge mark tail)
    resMatTmpCase3[resMat != 1] <- 0
    
    resMatTmpAll <- resMatTmpCase1 + resMatTmpCase2 + resMatTmpCase3
    resMatTmpAll[resMatTmpAll != 0] <- 1
    resMat <- resMatTmpAll
  }
  
  resMat
}

# pagIsAncestorEdge <- function(ancestralAmat, assumeNoSelectionVars = TRUE){
#   resMat <- ancestralAmat
#   # find edges of type
#   # 1)  amat[a,b] = 2  and  amat[b,a] = 3  implies   a --> b.
#   # 2) amat[a,b] = 1  and  amat[b,a] = 3   implies   a --o b.
#   
#   # Case 1: a --> b.
#   resMatTmpCase1 <- resMat + t(resMat)
#   resMatTmpCase1[resMatTmpCase1 != 5] <- 0
#   # remove entry so that resulting resMat[i,j] = 1 
#   # means that i is an ancestor of j
#   resMatTmpCase1[resMat != 2] <- 0
# 
#   if(assumeNoSelectionVars){
#     # Case 2: a --o b
#     resMatTmpCase2 <- resMat + t(resMat)
#     # retain edges containing circles
#     resMatTmpCase2[resMat != 1 & t(resMat) != 1] <- 0
#     # get only edges of type  a --o b
#     resMatTmpCase2[resMatTmpCase2 != 4] <- 0
#     # remove entry so that resulting resMat[i,j] = 1 
#     # means that i might be an ancestor of j
#     # -> only keep entries of type 1 so that
#     # amat[a,b] = 1 remain --> a might be an ancestor of b (b has edge mark circle,
#     # a has edge mark tail)
#     resMatTmpCase2[resMat != 1] <- 0
#     
#     # Case 3: a -- b (3,3)
#     resMatTmpCase3 <- resMat + t(resMat)
#     resMatTmpCase3[resMatTmpCase3 != 6] <- 0 
#     resMatTmpCase3[resMat != 3] <- 0
#   }else{
#     resMatTmpCase2 <- resMatTmpCase3 <- matrix(0, nrow = nrow(resMat), ncol = ncol(resMat))
#   }
#   
#   resMatTmpAll <- resMatTmpCase1 + resMatTmpCase2 + resMatTmpCase3
#   resMatTmpAll[resMatTmpAll != 0] <- 1
#   resMatTmpAll
#   
# }
# 
# pagIsMaybeAncestorEdge <- function(ancestralAmat, assumeNoSelectionVars = TRUE){
#   resMat <- ancestralAmat
#   # find edges of type
#   # amat[a,b] = 1  and  amat[b,a] = 1   implies   a o-o b.
#   # amat[a,b] = 1  and  amat[b,a] = 2   implies   a <-o b.
#   resMatTmp <- resMat + t(resMat)
#   # retain edges containing circles
#   resMatTmp[resMat != 1 & t(resMat) != 1] <- 0
#   
#   # Case 1: a o-o b
#   resMatTmpCase1 <- resMatTmp
#   # get only edges of type a o-o b
#   resMatTmpCase1[resMatTmpCase1 != 2] <- 0
#   
#   # Case 2: a o-> b
#   resMatTmpCase2 <- resMatTmp
#   # get only edges of type a <-o b 
#   resMatTmpCase2[resMatTmpCase2 != 3] <- 0
#   # remove entry so that resulting resMat[i,j] = 1 
#   # means that i might be an ancestor of j
#   # -> only keep entries of type 2 so that
#   # amat[b,a] = 2 remain --> b might be an ancestor of a (a has edge mark arrowhead,
#   # b has edge mark circle)
#   resMatTmpCase2[resMat != 2] <- 0
#   
#   #TODO add o--
#   
#   
#   if(!assumeNoSelectionVars){
#     # Case 3: a --o b
#     resMatTmpCase3 <- resMat + t(resMat)
#     # retain edges containing circles
#     resMatTmpCase3[resMat != 1 & t(resMat) != 1] <- 0
#     # get only edges of type  a --o b
#     resMatTmpCase3[resMatTmpCase3 != 4] <- 0
#     # remove entry so that resulting resMat[i,j] = 1 
#     # means that i might be an ancestor of j
#     # -> only keep entries of type 1 so that
#     # amat[a,b] = 1 remain --> a might be an ancestor of b (b has edge mark circle,
#     # a has edge mark tail)
#     resMatTmpCase3[resMat != 1] <- 0
#     
#   }else{
#     resMatTmpCase3 <- matrix(0, nrow = nrow(resMat), ncol = ncol(resMat))
#   }
#   
#   
#   resMatTmpAll <- resMatTmpCase1 + resMatTmpCase2 + resMatTmpCase3
#   resMatTmpAll[resMatTmpAll != 0] <- 1
#   resMatTmpAll
# }

pagIsMaybeAncestor <- function(ancestralAmat, assumeNoSelectionVars = TRUE){
  p <- ncol(ancestralAmat)
  cumsum <- matrix(0, p,p)
  paths <- lapply(1:p, function(j){ 
    bothPaths <- pagIsMaybeAncestorEdge(ancestralAmat, assumeNoSelectionVars) +
      pagIsAncestorEdge(ancestralAmat, assumeNoSelectionVars)
    
    bothPaths%^%j
  })
  for(i in 1:p) cumsum <- cumsum + paths[[i]]
  cumsum[cumsum != 0] <- 1
  diag(cumsum) <- 0
  cumsum
}


pagIsAncestor <- function(ancestralAmat, assumeNoSelectionVars = TRUE){
  p <- ncol(ancestralAmat)
  cumsum <- matrix(0, p,p)
  paths <- lapply(1:p, function(j) 
    pagIsAncestorEdge(ancestralAmat, assumeNoSelectionVars)%^%j)
  for(i in 1:p) cumsum <- cumsum + paths[[i]]
  cumsum[cumsum != 0] <- 1
  diag(cumsum) <- 0
  cumsum
}

pagIsNoAncestor <- function(ancestralAmat, assumeNoSelectionVars = TRUE){
  maybeAncestralMat <- pagIsMaybeAncestor(ancestralAmat, assumeNoSelectionVars)
  noAncestorMat <- 0*ancestralAmat
  # setting all entries to 1 that are neither ancestors nor possible ancestors
  noAncestorMat[maybeAncestralMat == 0] <- 1
  noAncestorMat
}
