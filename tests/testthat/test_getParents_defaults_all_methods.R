library(CompareCausalNetworks)
context("All supported methods")

data("simData_unknownShiftInterventions")

X <- simData_unknownShiftInterventions$X
environment <- simData_unknownShiftInterventions$environment

methods <- c("ICP", "hiddenICP", "backShift", "pc", "lingam", 
             "ges", "CAM", "rfci", "regression", 
             "bivariateANM", "bivariateCAM")


# TODO: modify title of check, put loop outside?
# TODO: change all method names to spelling in original package?
test_that("Checks configs for ICP", {
  for(method in methods){
    expect_is(
      Ahat <- getParents(X, environment, method=method, alpha=0.1)
      , "Matrix")
  }
})

# gies