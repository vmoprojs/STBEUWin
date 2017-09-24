####################################################
### Authors:  Moreno Bevilacqua, Víctor Morales Oñate
### Email:  moreno.bevilacqua@uv.cl, victor.morales@uv.cl
### Instituto de Estadistica
### Universidad de Valparaiso
### File name: STBEUTimeEvalOcl.r
### Description:
### This file contains a set of procedures
### for the set up of all the package routines.
### Last change: 13/09/2017.
####################################################
STBEUTimeEvalOcl<-function(theta,fix,coords,times,cc,data,type_dist,maxdist,maxtime, winc_s,winstp_s,
                   winc_t,winstp_t,subs,weighted,local=c(1,1),GPU=NULL)             
{
  model = cc
  names = checkpar(fix=fix,theta=theta,cc=model)
  path.parent <- getwd()
  if(!is.null(GPU)) 
  {
    if (model ==1) kernel = "DouExp.cl"
    if (model ==2) kernel = "Gneiting.cl"
    path <- system.file("CL", kernel, package = "STBEU")
    path <- gsub(paste("/",kernel,sep=""),"/",path);setwd(path)
    # fname <- paste(fname,"_OCL",sep="")
    .C("create_binary_kernel",  as.integer(GPU),as.character(kernel),  PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
  }
  
  setup=setting_param(cc,theta,fix);ntime=length(times);npts = nrow(coords)
  mom_cond = ww = grad = rep(0,setup$npar)
  gradcor = rep(0,setup$nparc);coordx=coords[,1];coordy=coords[,2];
  
  if(subs==1)
  {
    tCPU = proc.time()
    p = .C("scalar_space",as.integer(npts), as.integer(ntime),
           as.double(times), as.double(maxtime),as.double(maxdist),
           as.integer(cc), as.double(setup$parcor),as.integer(setup$flagcor),as.integer(setup$flagnuis),
           as.integer(setup$npar),as.double(setup$nuis),as.double(data),as.integer(weighted),
           mom_cond = as.double(mom_cond), as.integer(type_dist), as.double(coordx),
           as.double(coordy),as.double(gradcor),as.double(grad), ww= as.double(ww),  PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
    tCPU = proc.time()-tCPU
  }
  if(subs==2)
  {
    tCPU = proc.time()
    p = .C("scalar_time",as.integer(npts), as.integer(ntime),
           as.double(times), as.integer(cc),as.double(setup$parcor),as.integer(setup$flagcor),
           as.double(gradcor),as.integer(setup$flagnuis),as.double(grad),as.integer(setup$npar),
           as.double(setup$nuis),as.double(data),as.integer(weighted),as.double(maxtime),
           ww= as.double(ww), mom_cond = as.double(mom_cond), as.integer(type_dist),
           as.double(coordx),as.double(coordy), as.double(maxdist),PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
    tCPU = proc.time()-tCPU;tCPU
  }
  if(subs == 3)
  {
    tCPU = proc.time()
    p = .C("scalar_spacetime",as.integer(npts), as.integer(ntime),
           as.double(times), as.double(maxtime),
           as.integer(cc), as.double(setup$parcor),as.integer(setup$flagcor),
           as.double(gradcor),as.integer(setup$flagnuis),
           as.double(grad), as.integer(setup$npar),as.double(setup$nuis),as.double(data),as.integer(weighted), ww= as.double(ww),
           mom_cond = as.double(mom_cond),as.double(maxdist), as.integer(type_dist), as.double(coordx),
           as.double(coordy))
    tCPU = proc.time()-tCPU;tCPU
  }
  
  if(model==1){fname = "DoubleExpOCL"}
  if(model==2){fname = "GneitingOCL"}
  
  tGPU = proc.time()
  p1 =.C(fname,as.integer(npts), as.integer(ntime),
         as.double(times), as.double(maxtime),as.double(maxdist),
         as.integer(cc), as.double(setup$parcor),as.integer(setup$flagcor),as.integer(setup$flagnuis),
         as.integer(setup$npar),as.double(setup$nuis),as.double(data),as.integer(weighted),
         mom_cond = as.double(mom_cond), as.integer(type_dist), as.double(coordx),
         as.double(coordy),as.double(gradcor),as.double(grad), ww= as.double(ww),
         local = as.integer(local), GPU = as.integer(GPU),  PACKAGE='STBEU',DUP = TRUE, NAOK=TRUE)
  tGPU = proc.time()-tGPU 
  
  partime = list(ncoords = npts, ntime = ntime, maxdist = maxdist, maxtime = maxtime,
                 winc_s=winc_s,winstp_s=winstp_s,winc_t=winc_t,winstp_t=winstp_t,subs=subs) #time evaluation parameters
  if(tGPU[3]<tCPU[3]) ti = "OpenCL" else ti = "CPU"
  
  setwd(path.parent)
  return(list(cpupar_eval = p$mom_cond,openclpar_eval=p1$mom_cond,
         CPU_time=tCPU,OCL_time = tGPU, difftime = abs(tCPU-tGPU), faster = ti,partime = partime))
}