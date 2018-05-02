indtestAll <- function(fct,x,y,alpha,pars = list())
{
    result <- fct(x,y,alpha,pars)
}

indtestMutualAll <- function(f,x,alpha,pars = list())
{
    result <- f(x,alpha,pars)
}


indtestHsic <- function(x,y,alpha=0.05,pars)
{
    return(dhsic.test(x,y))
}
