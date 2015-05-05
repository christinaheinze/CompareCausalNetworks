
plotGraph <- function(A, main="",labels=NULL,layoutfunction=layout.circle,...){
    if(!is.matrix(A)) stop("A needs to be a matrix")
    if(nrow(A)!=ncol(A)) stop("A needs to have as many rows as columns")
    if(is.null(labels)) labels <- if( !is.null(cc <- colnames(A)))
        cc else as.character(1:ncol(A))
    G <- graph.adjacency(A,mode="directed",weighted="a")
    if(is.null(layoutfunction)) layoutfunction <-  layout.circle
    layout <- layoutfunction(G)

    optionals <- list(...)
    if(is.null(optionals$vertex.label.cex)) optionals$vertex.label.cex <- 1.5
    if(is.null(optionals$vertex.label.color)) optionals$vertex.label.color <-rgb(0.8,0.1,0.1,0.7)
    if(is.null(optionals$vertex.color)) optionals$vertex.color <-"white"
    if(is.null(optionals$vertex.frame.color)) optionals$vertex.frame.color <-rgb(0.8,0.1,0.1,0.5)
    if(is.null(optionals$edge.color)) optionals$edge.color <-rgb(0.1,0.1,0.1,0.5)
    if(is.null(optionals$edge.arrow.size)) optionals$edge.arrow.size <-0.7
    if(is.null(optionals$edge.arrow.width)) optionals$edge.arrow.width <-2
    if(is.null(optionals$vertex.size)) optionals$vertex.size <-30
    if(is.null(optionals$vertex.label.dist)) optionals$vertex.label.dist <-0
    if(is.null(optionals$vertex.label.degree)) optionals$vertex.label.degree <-  -pi/2

    plot(G, layout=layout, main=main, vertex.label=labels,vertex.shape="circle",vertex.label.cex=optionals$vertex.label.cex, vertex.label.color=optionals$vertex.label.color,vertex.color=optionals$vertex.color ,vertex.frame.color=optionals$vertex.frame.color , edge.color=optionals$edge.color, edge.arrow.size=optionals$edge.arrow.size ,edge.arrow.width=optionals$edge.arrow.width ,vertex.size=optionals$vertex.size,vertex.label.dist=optionals$vertex.label.dist,vertex.label.degree=optionals$vertex.label.degree)
    return(layout)
}
