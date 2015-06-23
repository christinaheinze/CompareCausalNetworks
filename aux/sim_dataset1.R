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
simData_unknownShiftInterventions <- simulateInterventions(n, p, A, G, intervMultiplier = 3, 
                                   noiseMult = 1, nonGauss = TRUE, hiddenVars = TRUE, 
                                   knownInterventions = FALSE, fracVarInt = NULL, simulateObs = TRUE, 
                                   seed = myseed)

save(simData_unknownShiftInterventions, file = "data/simData_unknownShiftInterventions.rda")
