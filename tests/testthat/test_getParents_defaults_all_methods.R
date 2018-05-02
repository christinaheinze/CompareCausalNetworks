library(CompareCausalNetworks)
library(testthat)
context("All supported methods")

data("simDataInv")

X <- simDataInv$X
environment <- simDataInv$environment
interventions <- simDataInv$interventions
mode <- "isAncestor"
methods <- c(
             "arges", "backShift", "bivariateANM", "bivariateCAM", "CAM",
             "fci", "fciplus", "ges", "gies", "hiddenICP", "ICP", "LINGAM",
             "mmhc", "rankArges", "rankFci", "rankGes", "rankGies", "rankPc",
             "rfci", "pc", "regression",
             "RESIT")


for(method in methods){
  test_that(paste("Checks output type for", method), {
    cat(paste("\nMethod:", method, "\n"))
    if(checkDependencies(method, checkRequireNamespace)){
      expect_is(
        Ahat <- getParents(X, environment, interventions, method=method, alpha=0.1, mode = mode, sparse = TRUE)
        , "Matrix")
      
      expect_is(
        Ahat <- getParents(X, environment, interventions, method=method, alpha=0.1, mode = mode, sparse = FALSE)
        , "matrix")
      
      print("\n")
      print(Ahat)
      
      expect_is(
        Ahat <- getParents(X, environment,interventions, method=method, alpha=0.1, mode = mode, returnAsList = TRUE)
        , "list")
    }else{
      print(paste("The required package for ", method," is not installed."))
    }
  }
  )
}