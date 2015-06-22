library(CompareCausalNetworks)
context("ICP")

require(backShift)

# 1st example:
# Simulate data with connectivity matrix A with assumptions 
# 1) hidden variables present
# 2) precise location of interventions is assumed unknown
# 3) different environments can be distinguished

## simulate data
myseed <- 1

# sample size n
n <- 1000

# p=3 predictor variables and connectivity matrix A
p  <- 3
A <- diag(p)*0
A[1,2] <- 0.8
A[2,3] <- 0.8
A[3,1] <- -0.4  

# divide data in 10 different environments
G <- 10

# simulate
simResult <- simulateInterventions(n, p, A, G, intervMultiplier = 3, 
noiseMult = 1, nonGauss = TRUE, hiddenVars = TRUE, 
knownInterventions = FALSE, fracVarInt = NULL, simulateObs = TRUE, 
seed = myseed)
X <- simResult$X
environment <- simResult$environment

## apply all  methods given in vector 'methods'
## (using all data pooled for pc/lingam/rfci -- can be changed with option 
## 'onlyObservationalData=TRUE')

method <- "rfci"

# Option 1): use this estimator as a point estimate
# Ahat <- getParents(X, environment, method=method, alpha=0.1)

# Option 2): use a stability selection based estimator
# with expected number of false positives bounded by EV=2
# Ahat <- getParentsStable(X, environment, EV=2, method=method, alpha=0.1)

test_that("Checks configs for ICP", {
  expect_is(
    Ahat <- getParents(X, environment, method=method, alpha=0.1)
    , "Matrix")
})


## simulate data
myseed <- 1

# sample size n
n <- 10000

# p=5 predictor variables and connectivity matrix A
p <- 5
A <- diag(p)*0
A[1,2] <- 0.8
A[2,3] <- -0.8
A[3,4] <- 0.8
A[3,5] <- 0.8
A[4,5] <- 0.3

# can add/remove feedback by using/not using
A[5,2] <- 0.8 

# divide data in 10 different environments
G <- 10

# simulate choose explicity intervention targets
simResult <- simulateInterventions(n, p, A, G, intervMultiplier = 3, 
 noiseMult = 1, nonGauss = TRUE, hiddenVars = FALSE, 
 knownInterventions = TRUE, fracVarInt = 0.1*p, simulateObs = TRUE, 
 seed = myseed)
X <- simResult$X
environment <- simResult$environment
interventions <- simResult$interventions

# number of unique environments
G <- length(unique(environment))

method <- "gies"

# Option 1): use this estimator as a point estimate
# Ahat <- getParents(X, environment, method=method, alpha=0.1)

# Option 2): use a stability selection based estimator
# with expected number of false positives bounded by EV=2
# Ahat <- getParentsStable(X, environment, EV=2, method=method, alpha=0.1)
# 
test_that("Checks configs for GIES", {
  expect_is(
    Ahat <- getParents(X, environment, interventions, method=method, alpha=0.1, 
                       excludeTargetInterventions = TRUE)
    , "Matrix")
})