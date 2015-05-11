\name{getParentsStable}
\alias{getParentsStable}
\title{
Estimate the connectivity matrix of a causal graph using stability selection
}
\description{
Estimates the connectivity matrix of a directed causal graph, using various possible methods. Supported methods at the moment are PC, Lingam, GES, GIES, RFCI, CAM, and invariant prediction with or without hidden variables (ICP and hiddenICP) as well as invariant covariance estimation (hiddenICE). Uses stability selection to select an appropriate sparseness.
}
\usage{
getParentsStable(X, environment, interventions = NULL, EV = 1, nodewise = TRUE, threshold = 0.75, 
                sampleSettings = 1/sqrt(2), sampleObservations = 1/sqrt(2), parentsOf = 1:ncol(X), 
                method = c("ICP", "hiddenICP", "hiddenICE", "pc", "lingam", "ges", "gies", "cam", "rfci")[1], 
                alpha = 0.1, variableSelMat = NULL, excludeTargetInterventions = TRUE, 
                onlyObservationalData = FALSE, indexObservationalData = 1,
                setOptions = list(), warnings = TRUE, nsim = 100)
}
\arguments{
  \item{X}{
A (nxp)-data matrix with n observations of p variables.
}
  \item{environment}{
An optional vector of length n, where the entry for observation i is an index for the environment in which observation i took place (simplest case entries \code{1} for observational data and entries \code{2} for interventional data of unspecified type). Is required for methods "ICP", "hiddenICP", "hiddenICE".}
  \item{interventions}{
A optional list of length n. The entry for observation i is a numeric vector that specifies the variables on which interventions happened for observation i (a scalar if an intervention happened on just one variable and \code{numeric(0)} if no intervention occured for this observation). Is used for method \code{gies} but will generate the vector \code{environment} if this is set to \code{NULL} (even though it might generate too many different environments for some data so a hand-picked vector \code{environment} is preferable). Is also used for \code{ICP} and \code{hiddenICP} to exclude interventions on the target variable of interest.}
  \item{EV}{
A bound on the expected number of falsely selected edges.
}
  \item{nodewise}{
Apply stability selection per node?
}
  \item{threshold}{
The empirical selection frequency in (0.5,1) under saubsampling that needs to be surpassed for an edge to be selected.
}
  \item{sampleSettings}{
The fraction of different environments to resample in each resampling (at least two different environments will be selected so the argument is without effect if there are just two different environments in total).
}
  \item{sampleObservations}{
The fraction of samples to resample in each environment.
}
  \item{parentsOf}{
The variables for which we would like to estimate the parents. Default are all variables.
}
  \item{method}{
A string that specfies the method to use. The methods \code{pc} (PC-algorithm), \code{lingam} (Lingam), \code{ges} (Greedy equivalence search), \code{gies} (Greedy interventional equivalence search) and \code{rfci} (Really fast causal inference) are imported from the package "pcalg" and are documented there in more detail, including the additional options that can be supplied via \code{setOptions}. The method \code{cam} (Causal additive models) is documented in the package "cam" and the methods \code{ICP} (Invariant causal prediction), \code{hiddenICP} (Invariant causal prediction with hidden variables) and \code{hiddenICE} (Invariant covariance estimation with hidden variables) are from the package "InvariantCausalPrediction". Finally, the methods \code{bivariateANM} and \code{bivariateCAM} are for now implemented internally but will hopefully be part of another package at some point in the near future.}
  \item{alpha}{
The level at which tests are done. This leads to confidence intervals for \code{ICP} and \code{hiddenICP} and is used internally for \code{pc} and \code{rfci}.}
  \item{variableSelMat}{
An optional logical matrix of dimension (pxp). An entry \code{TRUE} for entry (i,j) says that variable i should be considered as a potential parent for variable j and vice versa for \code{FALSE}. If the default value of \code{NULL} is used, all variables will be considered, but this can be very slow, especially for methods "pc", "ges", "gies", "rfci" and "cam".}
  \item{excludeTargetInterventions}{
When looking for parents of variable k in 1,...,p, set to \code{TRUE} if observations where an intervention on variable k occured should be excluded. Default is \code{TRUE}.
}
  \item{onlyObservationalData}{
If set to \code{TRUE}, only observational data is used. It will take the index in \code{environment} specified by \code{indexObservationalData}. Default is \code{FALSE}.
%If \code{interventions} is specified, it will read out observational data from this list, otherwise it will take the first occuring entry in \code{environment} as a code for observational data. Default is \code{FALSE}.
}
  \item{indexObservationalData}{
Index in \code{environment} that encodes observational data. Default is \code{1}.
}
  \item{setOptions}{
A list that can take method-specific options; see the individual documentations of the methods for more options and their possible values.
}
  \item{warnings}{
If set to \code{TRUE}, warnings will be printed. 
}
  \item{nsim}{
The number of resamples for stability selection.
}
}
\value{
A sparse matrix, where a 0 entry in (j,k) corresponds to an estimate of "no edge" \code{j} -> \code{parentsOf[k]}. Entries between 0 and 100 give the selection percentage of this edge over all resamples (set to 0 if below critical threshold) and all non-zero values are considered as selected edges.  
}
\references{
  Stability selection (2010):     
  N. Meinshausen and P. Buhlmann,
  Journal of the Royal Statistical Society: Series B, 72, 417-473
}
\author{
Nicolai Meinshausen <meinshausen@stat.math.ethz.ch>
}
\seealso{
\code{\link{getParents}} for the underlying point-estimate of the causal graph. \code{\link{plotGraph}} for plotting the estimated adjacency matrix if the whole graph is estimated. 
}
\examples{
}
\keyword{ Causality }
\keyword{ Graph estimation }
