####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate
### Email:  moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: SubSampTimeEval.r
### Description:
### This file contains a set of procedures
### for the set up of all the package routines.
### Last change: 13/09/2017.
####################################################
SubSampTimeEval<-function(theta,fix,coords,times,cc,data,type_dist,maxdist,maxtime, winc_s,winstp_s,
                   winc_t,winstp_t,subs,weighted,local=c(1,1),GPU=NULL)             
{
  names = checkpar(fix=fix,theta=theta,cc=cc)
  path.parent <- getwd()
  
  setup=setting_param(cc,theta,fix)
  # if(setup$parcor[1] <0 | setup$parcor[2]<0){
  #   obj = 1e100
  #   return (obj)}
  
  
  vari=double(setup$npar*0.5*(setup$npar-1)+setup$npar) ## vector of upper  triangular cov matrix (with diagonal)
  
  means=double(setup$npar)                               ## vector of means
  coordx=coords[,1];coordy=coords[,2];ncoords=nrow(coords);ntime=length(times)
  npar=length(theta)
  
  res1 = NULL;res2=NULL
  tGPU=NULL;tCPU=NULL
  if(!is.null(GPU)) 
  {
    if (model ==1) kernel = "DouExp.cl"
    if (model ==2) kernel = "Gneiting.cl"
    path <- system.file("CL", kernel, package = "STBEU")
    path <- gsub(paste("/",kernel,sep=""),"/",path);setwd(path)
    # fname <- paste(fname,"_OCL",sep="")
    .C("create_binary_kernel",  as.integer(GPU),as.character(kernel),  PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
    
    if(subs==1 && cc==1) type_sub="SubSamp_space_DE_ocl"
    if(subs==2 && cc==1) type_sub="SubSamp_time_DE_ocl"
    if(subs==3 && cc==1) type_sub="SubSamp_spacetime_DE_ocl"
    
    if(subs==1 && cc==2) type_sub="SubSamp_space_GN_ocl"
    if(subs==2 && cc==2) type_sub="SubSamp_time_GN_ocl"
    if(subs==3 && cc==2) type_sub="SubSamp_spacetime_GN_ocl"
    
    
    tGPU = proc.time()
    res1=.C(type_sub,as.double(coordx), as.double(coordy), as.double(times),as.integer(ncoords),as.integer(ntime),as.integer(cc), as.double(data),as.integer(type_dist),as.double(maxdist),as.double(maxtime),as.integer(setup$npar),as.double(setup$parcor),as.integer(setup$nparc),
         as.double(setup$nuis),as.integer(setup$nparnuis),
         as.integer(setup$flagcor),as.integer(setup$flagnuis),vv=as.double(vari),as.double(winc_s), as.double(winstp_s),
         as.double(winc_t), as.double(winstp_t),
         mm=as.double(means),as.integer(weighted),as.integer(local),as.integer(GPU))
    tGPU = proc.time()-tGPU;tGPU
    
    # x=res1$mm       #means vector
    # F=xpnd(res1$vv) #cov matrix
    # Fchol = MatDecomp(F,"cholesky")
    # if(length(Fchol)==1)
    # {
    #   Fchol = MatDecomp(F,"svd")
    #   inv = MatInv(Fchol,"svd")
    #   obj1= crossprod(crossprod(inv, x), x)
    # }
    # else
    # {
    #   inv=chol2inv(Fchol)  #inverse   
    #   obj1= crossprod(crossprod(inv, x), x)# quadratic form
    # }
  }
  {
    if(subs==1) type_sub="SubSamp_space"
    if(subs==2) type_sub="SubSamp_time"
    if(subs==3) type_sub="SubSamp_spacetime"
    
    tCPU = proc.time()
    res2=.C(type_sub,as.double(coordx), as.double(coordy), as.double(times),as.integer(ncoords),as.integer(ntime),as.integer(cc), as.double(data),as.integer(type_dist),as.double(maxdist),as.double(maxtime),as.integer(setup$npar),as.double(setup$parcor),as.integer(setup$nparc),
         as.double(setup$nuis),as.integer(setup$nparnuis),
         as.integer(setup$flagcor),as.integer(setup$flagnuis),vv=as.double(vari),as.double(winc_s), as.double(winstp_s),
         as.double(winc_t), as.double(winstp_t),
         mm=as.double(means),as.integer(weighted),as.integer(local),as.integer(GPU))
    tCPU = proc.time()-tCPU;tCPU
    
    # x=res2$mm       #means vector
    # #print(x)
    # F=xpnd(res2$vv) #cov matrix
    # #print(F)
    # #F <- matrix(c(-1,0,0,1),ncol=2)
    # Fchol = MatDecomp(F,"cholesky")
    # if(length(Fchol)==1)
    # {
    #   #print(F)
    #   Fchol = MatDecomp(F,"svd")
    #   inv = MatInv(Fchol,"svd")
    #   obj2= crossprod(crossprod(inv, x), x)
    # }
    # else
    # {
    #   
    #   inv=chol2inv(Fchol)  #inverse   
    #   obj2= crossprod(crossprod(inv, x), x)# quadratic form
    # }
  }
 setwd(path.parent)
 partime = list(ncoords = ncoords, ntime = ntime, maxdist = maxdist, maxtime = maxtime,
                winc_s=winc_s,winstp_s=winstp_s,winc_t=winc_t,winstp_t=winstp_t,subs=subs) #time evaluation parameters
 if(tGPU[3]<tCPU[3]) ti = "OpenCL" else ti = "CPU"
 return(list(cpupar_eval = res2$mm,openclpar_eval=res1$mm,
             CPU_time=tCPU,OCL_time = tGPU, difftime = abs(tCPU-tGPU), faster = ti,partime=partime))
}