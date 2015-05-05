\name{plotGraph}
\alias{plotGraph}
\title{
A basic plotting function to visualize directed graphs
}
\description{
Given an estimated adjacency matrix, the function plots the corresponding graph.
}
\usage{
plotGraph(A, main = "", labels = NULL, layoutfunction = layout.circle,...)
}
\arguments{
  \item{A}{
A (pxp) adjacency matrix.
}
  \item{main}{
A title for the graph
}
  \item{labels}{
Names for the variabels. If left at the default \code{NULL}, it uses the column-names from the adjacency matrix if they are non-null.
}
  \item{layoutfunction}{
The layout-function used for plotting the graph. Defaults to \code{layout.circle} but other options (see igraph-package) can be used such as \code{layout.kamada.kawai} etc.
}
  \item{...}{
Optional arguments passed to the plotting function. Consists of igraph-type options like \code{vertex.label.cex},\code{vertex.label.color}, \code{edge.arrow.size} or \code{vertex.size} etc.
}
}
\details{
Currently not all options of igraph are used; additional arguments are ignored. 
}
\value{
Returns the layout of the graph which can be used for further plots.
}

\author{
Nicolai Meinshausen <meinshausen@stat.math.ethz.ch>
}
\seealso{
\code{\link{getParents}} and \code{\link{getParentsStable}} for the estimation of the causal graph.
}