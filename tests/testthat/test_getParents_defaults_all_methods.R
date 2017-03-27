library(CompareCausalNetworks)
context("All supported methods")

data("simDataInv")

X <- simDataInv$X
environment <- simDataInv$environment
interventions <- simDataInv$interventions
mode <- "isAncestor"
methods <- c("arges", "backShift", "bivariateANM", 
             "bivariateCAM", "CAM", 
             "fci", "fciplus", "ges", "gies", "hiddenICP",
             "ICP", "LINGAM", "mmhc", "rankArges", "rankFci",
             "rankGes", "rankGies", "rankPc", "rfci", "pc",
             "regression")


# TODO: change all method names to spelling in original package?

for(method in methods){
  test_that(paste("Checks output type for", method), {
    
    expect_is(
      Ahat <- getParents(X, environment, interventions, method=method, alpha=0.1, mode = mode, sparse = TRUE)
      , "Matrix")
    
    expect_is(
      Ahat <- getParents(X, environment, interventions, method=method, alpha=0.1, mode = mode, sparse = FALSE)
      , "matrix")
    
    cat(paste("\nMethod:", method, "\n"))
    print(Ahat)
    
    expect_is(
      Ahat <- getParents(X, environment,interventions, method=method, alpha=0.1, mode = mode, returnAsList = TRUE)
      , "list")
    
  }
  )
}