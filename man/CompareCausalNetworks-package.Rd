\name{CompareCausalNetworks-package}
\alias{CompareCausalNetworks-package}
\alias{CompareCausalNetworks}
\docType{package}
\title{
Compare estimates of causal graphs using a unified interface to various methods
}
\description{
Provides a unified interface to various causal graph estimation methods.
}
\details{
\tabular{ll}{
Package: \tab CompareCausalNetworks\cr
Type: \tab Package\cr
Version: \tab 0.2.2\cr
Date: \tab 2018-05-18\cr
License: \tab GPL\cr
}
The causal graphs can be estimated with function \code{\link{getParents}} 
and a stability-selection version is available at \code{\link{getParentsStable}}. 

The supported methods are provided through the packages listed in \code{Suggests}. 
Thus, to use a particular method the corresponding package needs to be installed 
on your machine. To run the examples, most of these packages need to be installed.
}
\author{
Christina Heinze-Deml <heinzedeml@stat.math.ethz.ch>, 
Nicolai Meinshausen <meinshausen@stat.math.ethz.ch>
}