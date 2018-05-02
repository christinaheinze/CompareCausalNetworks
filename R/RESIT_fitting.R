# Copyright (c) 2010-2013  Jonas Peters  [peters@stat.math.ethz.ch]
# All rights reserved.  See the file COPYING for license terms. 
# 
# train_*** <- function(X, y, pars)
# e.g. train_linear or train_gp
#
# Performs a regression from X to y
# 
# INPUT:
#   X         nxp matrix of training inputs (n data points, p dimensions)
#   Y         vector of N training outputs (n data points)
#   pars      list containing parameters of the regression method
#
# OUTPUT:
#   result    list with the result of the regression
#      $model        list of learned model (e.g., weight vector)
#      $Yfit         fitted outputs for training inputs according to the learned model
#      $residuals    noise values (e.g., residuals in the additive noise case)


####
#Linear Regression
####
train_linear <- function(X,y,pars = list())
{
    mod <- lm(y ~ X)
    result <- list()
    result$Yfit = as.matrix(mod$fitted.values)
    result$residuals = as.matrix(mod$residuals)
    result$model = mod
    #for coefficients see list(mod$coef)
    return(result)
}

####
#gam Regression from mgcv package
####
train_gam <- function(X,y,pars = list(numBasisFcts = 10))
{
    if(!("numBasisFcts" %in% names(pars) ))
    { 
        pars$numBasisFcts = 10
    }
    p <- dim(as.matrix(X))
    if(p[1]/p[2] < 3*pars$numBasisFcts)
    {
        pars$numBasisFcts <- ceiling(p[1]/(3*p[2]))
        cat("changed number of basis functions to    ", pars$numBasisFcts, "    in order to have enough samples per basis function\n")
    }
    dat <- data.frame(as.matrix(y),as.matrix(X))
    coln <- rep("null",p[2]+1)
    for(i in 1:(p[2]+1))
    {
        coln[i] <- paste("var",i,sep="")
    }
    colnames(dat) <- coln
    labs<-"var1 ~ "
    if(p[2] > 1)
    {
        for(i in 2:p[2])
        {
            labs<-paste(labs,"s(var",i,",k = ",pars$numBasisFcts,") + ",sep="")
            #      labs<-paste(labs,"s(var",i,") + ",sep="")
            # labs<-paste(labs,"lo(var",i,") + ",sep="")
        }
    }
    labs<-paste(labs,"s(var",p[2]+1,",k = ",pars$numBasisFcts,")",sep="")
    # labs<-paste(labs,"s(var",p[2]+1,")",sep="")
    # labs<-paste(labs,"s(var",p[2]+1,", bs = "cc")",sep="") #factor 2 faster
    # labs<-paste(labs,"s(var",p[2]+1,", bs = "cr")",sep="") # factor 2 + eps faster
    # labs<-paste(labs,"lo(var",p[2]+1,")",sep="")
    mod_gam <- FALSE
    try(mod_gam <- gam(formula=formula(labs), data=dat),silent = TRUE)
    if(typeof(mod_gam) == "logical")
    {
        cat("There was some error with gam. The smoothing parameter is set to zero.\n")
        labs<-"var1 ~ "
        if(p[2] > 1)
        {
            for(i in 2:p[2])
            {
                labs<-paste(labs,"s(var",i,",k = ",pars$numBasisFcts,",sp=0) + ",sep="")
            }
        }
        labs<-paste(labs,"s(var",p[2]+1,",k = ",pars$numBasisFcts,",sp=0)",sep="")
        mod_gam <- gam(formula=formula(labs), data=dat)
    }
    result <- list()
    result$Yfit <- as.matrix(mod_gam$fitted.values)
    result$residuals <- as.matrix(mod_gam$residuals)
    result$model <- mod_gam 
    result$df <- mod_gam$df.residual     
    result$edf <- mod_gam$edf     
    result$edf1 <- mod_gam$edf1     
    
    # for degree of freedom see mod_gam$df.residual
    # for aic see mod_gam$aic
    return(result)
}

####
#gam boost regression from package mboost
####
train_GAMboost <- function(X,y,pars = list()) #
{
    ## begin old version
    # EXPLANATION: surprisingly, it turned out that this cannot be applied to large p (private discussion with T. Hothorn in Sep 2013)
    # yy <- y
    # dat <- data.frame(cbind(yy,X))
    # gb <- gamboost(yy ~ .,data=dat, baselearner = "bbs")
    ## end old version
    
    ## begin new version
    dat <- as.data.frame(X)
    bl <- lapply(dat, bbs)
    gb <- mboost_fit(bl, y)
    ## end new version
    
    result <- list()
    result$Yfit <- gb$fitted()
    result$residuals <- gb$resid()
    result$model <- gb
    return(result)
}

####
#lasso regression from package glmnet
####
train_lasso <- function(X,y,pars = list())
{
    cvres <- cv.glmnet(X,y)
    mod <- glmnet(X,y,lambda = cvres$lambda.1se)
    result <- list()
    result$Yfit <- predict(mod,X)
    result$residuals <- y - result$Yfit
    result$model <- mod
    return(result)
}

# =========
# 2. train_model
# =========

train_model <- function(f,X,y,pars = list())
{
    result <- f(X,y,pars)
}