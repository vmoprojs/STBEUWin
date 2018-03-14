####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate
### Email:  moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: Utility.r
### Description:
### This file contains a set of procedures
### for the set up of all the package routines.
### Last change: 13/09/2017.
####################################################

### decomposition of a square  matrix
MatDecomp<-function(mtx,method)    {
  # print(mtx)
  if(method=="cholesky")  {
    mat.decomp <- try(chol(mtx), silent=TRUE)
    if (inherits(mat.decomp , "try-error")) return (FALSE)
  }
  if(method=="svd")      {
    mat.decomp <- svd(mtx, nv = 0)
    cov.logdeth <- try(sum(log(sqrt(mat.decomp$d))), silent=TRUE)
    if (inherits(cov.logdeth, "try-error"))  return (FALSE)
  }
  return(mat.decomp)
}
### inverse a square matrix given a decomposition
MatInv<-function(mat.decomp,method)    {
  
  if(method=="cholesky")  varcov.inv <- chol2inv(mat.decomp)
  if(method=="svd")       
  { 
    tol = sqrt(.Machine$double.eps)
    e <- mat.decomp$d;e[e > tol] <- 1/e[e > tol] 
    varcov.inv<-mat.decomp$u %*% diag(e,nrow=length(e)) %*% t(mat.decomp$u) 
  } 
  return(varcov.inv)
} 

####################################################################################################################
####################################################################################################################


DevOpenCL <- function()
{
    .C("DevOpenCL",PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
}

#######################################################################
setting_param<-function(cc,theta,fix)
{
  param <- c(theta,fix)
  if(cc==1){      ### double exponential model
    npar=length(theta)
    parcor=c(param["scale_s"],param["scale_t"])      #corr parameters
    flagcor=c(1,1) 
    nuis=c(param["mean"],param["nugget"],param["sill"])   #mean, nugget, sill
    #flagnuis=c(1,0,1)
    flagnuis=c(0,0,1)
    nparc=sum(flagcor)  
    nparnuis=sum(flagnuis)
  }
  if(cc==2){       ### gneiting model
    npar=length(theta)
    parcor=c(param["power_s"],param["power_t"],param["scale_s"],param["scale_t"],param["sep"])     #corr parameters
    flagcor=c(0,0,1,1,0) 
    nuis=c(param["mean"],param["nugget"],param["sill"])   #mean, nugget, sill
    # flagnuis=c(1,0,1)
    flagnuis=c(0,0,1) 
    nparc=sum(flagcor)  
    nparnuis=sum(flagnuis)
  }
  setup=list(flagcor=flagcor,flagnuis=flagnuis,npar=npar,parcor=parcor,nuis=nuis,nparc=nparc
             ,nparnuis=nparnuis)
  return(setup)
}



checkpar <- function(fix,theta,cc)
{
  if(cc==1)
  {
    namespar = c("nugget","mean","scale_t","scale_s","sill")
    if((length(fix)+length(theta))!=5) stop("Number of parameters don't match the declared model")
    if(is.null(names(fix)) || is.null(names(theta))) 
    {
      warning("fix and theta should be named. If not, take care about the ouput order") 
    }else
    {
      namesfix = names(fix);namestheta = names(theta)
    }
    res = names(theta)
  }
  if(cc==2)
  {
    namespar = c("nugget","mean","scale_t","scale_s","sill",
                 "power_s","power_t","sep")
    if((length(fix)+length(theta))!=8) stop("Number of parameters don't match the declared model")
    if(is.null(names(fix)) || is.null(names(theta))) 
    {
      warning("fix and theta should be named. If not, take care about the ouput order") 
    }else
    {
      namesfix = names(fix);namestheta = names(theta)
    }
    res = names(theta)
  }
  return(res)
}
# checkpar(fix,theta)




print.STBEUFit <- function(x,names,GPU, digits = max(3, getOption("digits") - 3), ...)
{
  cat('\n=================================================================')
  if(x$convergence==0) cat('\nSuccessful convergence') else cat('\nConvergence may have failed')
  if(is.null(GPU)) cat('\nCPU device for computation') else cat('\nOpenCL device for computation')
  if(is.null(names)) {cat('\n Parameters (order of theta)\n',x$par)} else
  {
    res = x$par;names(res) = names
    cat('\nParameters:\n')
    print.default(res)
  }
  cat('\n=================================================================\n')
}
