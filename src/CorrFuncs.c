#include "header.h"


//#define EPS 1.0e-60

//---------START WENDLAND FUNCTIONS-----------

/* START Wendland covariance */

/* integrand  in  generalized wendland function*/
double int_gen(double x,double mu, double alpha,double lag,double supp)
{
    double res=0.0,y;
    y=lag/supp;
    res=R_pow(1-x,mu-1)*R_pow(x*x-y*y,alpha)/beta(2*alpha+1,mu);
    return (res);///(R_pow(2,alpha-1)*gamma(alpha)*R_pow(supp,2*alpha)));
}
void integr_gen(double *x, int n, void *ex){
    int i;double mu,alpha,beta,y;
    mu =    ((double*)ex)[0];  //mu
    alpha = ((double*)ex)[1];  //alpha
    beta =     ((double*)ex)[2];  //csupp
    y =     ((double*)ex)[3];  //h
    for (i=0;i<n;i++) {x[i]=int_gen(x[i],mu,alpha,y,beta);}
    return;
}
// function computing generalized wendland
double wendintegral(double x, double *param) {
    
    double ex[4], lower, upper, epsabs, epsrel, result, abserr, *work;
    int neval, ier, subdiv, lenw, last, *iwork;
    subdiv = 100;
    epsabs = R_pow(DOUBLE_EPS, 0.25);
    epsrel = epsabs;
    lenw = 4 * subdiv;           /* as instructed in WRE */
    iwork =   (int *) Calloc(subdiv, int);  /* idem */
    work = (double *) Calloc(lenw, double); /* idem */
    ex[0] = param[0]; ex[1] = param[1]; ex[2] = param[2];ex[3]=x;
    lower=x/param[2];
    upper=1;
    // Compute the integral
    
    if(x<=param[2]) {
        Rdqags(integr_gen, (void *) &ex,
               &lower, &upper, &epsabs, &epsrel, &result,
               &abserr, &neval, &ier, &subdiv, &lenw, &last, iwork, work);
        
    }else   {result=0;}
    Free(iwork);Free(work);
    return(result);
}
/* generalized wendland function*/
double CorFunW_gen(double lag,double R_power1,double smooth,double scale)  // mu alpha beta
{
    
    double rho=0.0,x=0;
    if(smooth==0) {
        
        x=lag/scale;
        if(x<=1) rho=R_pow(1-x,R_power1);
        else rho=0;
        return(rho);
    }
    if(smooth==1) {
        
        x=lag/scale;
        if(x<=1) rho=R_pow(1-x,R_power1+1)*(1+x*(R_power1+1));
        else rho=0;
        return(rho);
    }
    if(smooth==2) {
        
        x=lag/scale;
        if(x<=1) rho=R_pow(1-x,R_power1+2)*(1+x*(R_power1+2)+x*x*(R_power1*R_power1 +4*R_power1 +3 )/3  );
        else rho=0;
        return(rho);
    }
    
    x=lag;
    double *param;
    param=(double *) Calloc(3,double);
    param[0]=R_power1;param[1]=smooth;param[2]=scale;  //mu,alpha //beta
    
    rho=wendintegral(x,param);
    Free(param);
    return(rho);
}


double wen_time(double *par, double h,double u)
{
    double R_power_s=2.0;
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];

    double arg=R_pow(1+R_pow(h/scale_s,R_power_s/2),-1/(R_power_s/2));
    double rho=R_pow(arg,R_power)*CorFunW_gen(u,R_power_t,smooth,scale_t*R_pow(arg,sep));
    //printf("\n**** rho: %f\n",rho);
    return(rho);
    
    //return(0.5);
}



/* END Wendland covariance */


/* START DERIVATIVES Wendland covariance */

// SCALE_S:
double deri_scale_s_wen_time(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*scale_s;
    double delta=SQE*scale_s;
    
    double *par1;
    par1 = (double *) Calloc(7, double);
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s + delta;
    par1[4]=scale_t;
    par1[5]=sep;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //printf("+++ deri_scale_s_wen_time   %f\n",grad);
    //printf("+++ elementWise:   %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta);
     //printf("+++ scale_s   %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta);
    Free(par1);
    
    /*double paro11[7];
    paro11[0]=par[0];
    paro11[1]=par[1];
    paro11[2]=par[2];
    paro11[3]=par[3];
    paro11[4]=par[4];
    paro11[5]=par[5];
    paro11[6]=par[6];

    double delta=SQE*par[3];
    
    double paro1[7];
    paro1[0]=par[0];
    paro1[1]=2;
    paro1[2]=par[2];
    paro1[3]=par[3] + delta;
    paro1[4]=par[4];
    paro1[5]=par[5];
    paro1[6]=par[6];
    double grad=(wen_time(paro1[0],paro1[1],paro1[2],paro1[3],paro1[4],paro1[5],paro1[6],h,u)-wen_time(paro11[0],paro11[1],paro11[2],paro11[3],paro11[4],paro11[5],paro11[6],h,u))/delta;*/
    
    //printf("+++ deri_scale_s_wen_time   %f\n",grad);
    
   

    return(grad);
}



// SCALE_T:
double deri_scale_t_wen_time(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*scale_t;
    double delta=SQE*scale_t;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s;
    par1[4]=scale_t+ delta;
    par1[5]=sep;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //printf("+++ deri_scale_t_wen_time   %f\n",grad);
    //printf("+++ elementWise:   %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta);
    
    //printf("+++ scale_t   %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta);
    Free(par1);
    return(grad);
    
}



// SMOOTH:
double deri_smooth_wen_time(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*smooth;
    double delta=SQE*smooth;
    
    double *par1;
    par1 = (double *) Calloc(7, double);
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s;
    par1[4]=scale_t;
    par1[5]=sep;
    par1[6]=smooth+ delta;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //printf("+++ deri_smooth_wen_time   %f\n",grad);
    //printf("+++ elementWise:   %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta);
   Free(par1);
    return(grad);
}



// SILL:
double deri_sill_wen_time(double *par, double h,double u)
{
    double grad=wen_time(par,h,u);
    return(grad);
}



// SEP:
double deri_sep_wen_time(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*sep;
    double delta=SQE*sep;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s;
    par1[4]=scale_t;
    par1[5]=sep+ delta;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //printf("+++ sep   %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta);
    //printf("+++ deri_sep_wen_time   %f\n",grad);
   Free(par1);
    return(grad);
}



// R_power:
double deri_R_power_wen_time(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*R_power;
    double delta=SQE*R_power;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power+ delta;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s;
    par1[4]=scale_t;
    par1[5]=sep;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //printf("+++ deri_R_power_wen_time   %f\n",grad);
   Free(par1);
    return(grad);
}



// R_power_t:
double deri_R_power_t_wen_time(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*R_power_t;
    double delta=SQE*R_power_t;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t+ delta;
    par1[3]=scale_s;
    par1[4]=scale_t;
    par1[5]=sep;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
   // printf("+++ deri_R_power_t_wen_time   %f\n",grad);
    Free(par1);
    return(grad);
}




// NUGGET:
/*
double deri_nugget_wen_time(double *par, double h,double u)
{
    double mean=par[0],nugget=par[1],sill=par[2];
    
    (wen_time(par1,h,u)-wen_time(par,h,u))/delta
    double grad=wen_time(par,h,u);
    return(grad);
}
 */

/* END DERIVATIVES Wendland covariance */


/* START DERIVATIVES Wendland covariance */

// SCALE_S:
double deri_scale_s_wen_time_v(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*scale_s;
    double delta=SQE1*scale_s;
    
    double *par1;
    par1 = (double *) Calloc(7, double);
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s + delta;
    par1[4]=scale_t;
    par1[5]=sep;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //printf("+++ deri_scale_s_wen_time   %f\n",grad);
    //printf("+++ elementWise:   %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta);
     //printf("+++ scale_s   %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta);
    Free(par1);
    
    /*double paro11[7];
    paro11[0]=par[0];
    paro11[1]=par[1];
    paro11[2]=par[2];
    paro11[3]=par[3];
    paro11[4]=par[4];
    paro11[5]=par[5];
    paro11[6]=par[6];

    double delta=SQE*par[3];
    
    double paro1[7];
    paro1[0]=par[0];
    paro1[1]=2;
    paro1[2]=par[2];
    paro1[3]=par[3] + delta;
    paro1[4]=par[4];
    paro1[5]=par[5];
    paro1[6]=par[6];
    double grad=(wen_time(paro1[0],paro1[1],paro1[2],paro1[3],paro1[4],paro1[5],paro1[6],h,u)-wen_time(paro11[0],paro11[1],paro11[2],paro11[3],paro11[4],paro11[5],paro11[6],h,u))/delta;*/
    
    //printf("+++ deri_scale_s_wen_time   %f\n",grad);
    
   

    return(grad);
}



// SCALE_T:
double deri_scale_t_wen_time_v(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*scale_t;
    double delta=SQE1*scale_t;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s;
    par1[4]=scale_t+ delta;
    par1[5]=sep;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //printf("+++ deri_scale_t_wen_time   %f\n",grad);
    //printf("+++ elementWise:   %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta);
    
    //printf("+++ scale_t   %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta);
    Free(par1);
    return(grad);
    
}



// SMOOTH:
double deri_smooth_wen_time_v(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*smooth;
    double delta=SQE1*smooth;
    
    double *par1;
    par1 = (double *) Calloc(7, double);
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s;
    par1[4]=scale_t;
    par1[5]=sep;
    par1[6]=smooth+ delta;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //printf("+++ deri_smooth_wen_time   %f\n",grad);
    //printf("+++ elementWise:   %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta);
   Free(par1);
    return(grad);
}



// SILL:
double deri_sill_wen_time_v(double *par, double h,double u)
{
    double grad=wen_time(par,h,u);
    return(grad);
}



// SEP:
double deri_sep_wen_time_v(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*sep;
    double delta=SQE1*sep;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s;
    par1[4]=scale_t;
    par1[5]=sep+ delta;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //printf("+++ sep   %f %f %f\n",wen_time(par1,h,u),wen_time(par,h,u),delta);
    //printf("+++ deri_sep_wen_time   %f\n",grad);
   Free(par1);
    return(grad);
}



// R_power:
double deri_R_power_wen_time_v(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*R_power;
    double delta=SQE1*R_power;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power+ delta;
    par1[1]=2;
    par1[2]=R_power_t;
    par1[3]=scale_s;
    par1[4]=scale_t;
    par1[5]=sep;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
    //printf("+++ deri_R_power_wen_time   %f\n",grad);
   Free(par1);
    return(grad);
}



// R_power_t:
double deri_R_power_t_wen_time_v(double *par, double h,double u)
{
    double R_power=par[0];
    double R_power_t=par[2];
    double scale_s=par[3];
    double scale_t=par[4];
    double sep=par[5];
    double smooth=par[6];
    
    //double EPS = 1.0e-10;
    
    //double delta=sqrt(EPS1)*R_power_t;
    double delta=SQE1*R_power_t;
    
    double *par1;
    par1 = (double *) Calloc(7, double); /* idem */
    // define vector par1
    par1[0]=R_power;
    par1[1]=2;
    par1[2]=R_power_t+ delta;
    par1[3]=scale_s;
    par1[4]=scale_t;
    par1[5]=sep;
    par1[6]=smooth;
    double grad=(wen_time(par1,h,u)-wen_time(par,h,u))/delta;
   // printf("+++ deri_R_power_t_wen_time   %f\n",grad);
    Free(par1);
    return(grad);
}





//---------END WENDLAND FUNCTIONS-----------







double log_biv_Norm(double corr,double zi,double zj,double mi,double mj,double vari, double nugget)
{
    double u,v,u2,v2,det,s1,s12,dens;
    u=zi-mi;
    v=zj-mj;
    u2=R_pow(u,2);v2=R_pow(v,2);
    s1=vari;
    s12=vari*corr;
    det=R_pow(s1,2)-R_pow(s12,2);
    dens=(-0.5*(2*log(2*M_PI)+log(det)+(s1*(u2+v2)-2*s12*u*v)/det));
return(dens);
}




/********************************************************************/
/************ gradient of the pairwise CL ****************/
/********************************************************************/
void Grad_Pair_Gauss(double rho, int *flag, double *gradcor, double *grad,
                     int *npar, double *par, double u, double v,
                     int *cormod,double lags,double lagt,int *flagcor, int *nparcT,double *parcor)
{
    // Initialization variables:
    int j=0;
    double mean=par[0],nugget=par[1],sill=par[2], *parC;
    
    parC=(double *) Calloc(nparcT[0],double);
    //printf("------- npar[0]: %d nparcT[0]: %d\n",npar[0],nparcT[0]); 4, 6
    for(j=0;j<nparcT[0];j++) parC[j]=parcor[j];
   
    double a=sill,b=sill*rho*(1-nugget),pa=a*a,pb=b*b;
    double c=-pa+pb,d=pa+pb,k=1/(c*c);
    double C=0.0,L=0.0,R=0.0, delta=0.0,rhod,zi,zj;
    double pu=0.0, pv=0.0;
    int h=0, i=0,  kk;
    //defines useful quantities:
    zi=u-mean;
    zj=v-mean;
    pu=pow(zi,2);
    pv=pow(zj,2);
    R=pu+pv;L=zi*zj;
    
    double ff=log_biv_Norm(rho*(1-nugget),u,v, mean, mean,sill,nugget);
    // Derivatives  respect with the mean
    if(flag[0]==1)
    {
        grad[i]= (((zi*a - zj*b)) + ((zj*a - zi*b)))/(c);i++;
        
        
    }
    // Derivative  respect with the nugget
    if(flag[1]==1)
    {
    delta=SQE*nugget;
    grad[i]=(log_biv_Norm(rho*(1-(nugget+delta)),u,v, mean, mean,sill,0) -ff)/(delta);
        i++;
        //printf("nugget grad: %f\n",grad[i]);
       // grad[i]=0.5*k*(R*d-L*4*b*a-2*a*(pa-pb));i++;
    }
    // Derivative respect with the sill
    if(flag[2]==1)
    {
        delta=SQE*sill;
        grad[i]=(log_biv_Norm(rho*(1-(nugget)),u,v, mean, mean,sill+delta,0) -ff)/(delta);
        //printf("sill grad: %f log: %f ff: %f SQE: %f nugget %f delta: %f\n",grad[i],log_biv_Norm(rho,u,v, mean, mean,sill+delta,nugget),ff,SQE, nugget,delta);
        //printf("sill grad: %f\n",grad[i]);
        i++;
    }
    // Derivatives with respect to the correlation parameters
    
    /*
    h=0;
    C=-k*sill*(R*a*b-L*d+b*c);
    for(j=i;j<*npar;j++){
        //printf("cr grad: %f \n",grad[i]);
        grad[j]=C*gradcor[h];h++;
        
    }*/
     
    h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) {
     
  if(flagcor[h]==1){
       delta=SQE*parcor[h];
       parC[h]=parcor[h]+delta;
       rhod=CorFct(cormod,lags,lagt,parC);
      
      //printf("kk+i: %d h: %d rhod: %f lags: %f lagt: %f parC[h] %f parC[0] %f parC[1] %f parC[2] %f parC[3] %f parC[4] %f parC[5] %f parC[6] %f \n",kk+i,h,rhod,lags,lagt,parC[h],parC[0],parC[1],parC[2],parC[3],parC[4],parC[5],parC[6]);
      
     /* if(ISNAN(rhod))
      {
          printf("h: %d rhod: %f lags: %f lagt: %f parC[h] %f parC[0] %f parC[1] %f parC[2] %f parC[3] %f parC[4] %f parC[5] %f parC[6] %f \n",h,rhod,lags,lagt,parC[h],parC[0],parC[1],parC[2],parC[3],parC[4],parC[5],parC[6]);
      }*/
       
    grad[kk+i]=(log_biv_Norm(rhod*(1-nugget),u,v, mean, mean,sill,0) - ff)/delta;
      //printf("h: %d delta: %f parC[h]: %f rhod %f ff %f grad[kk+i]: %f\n",h, delta,parC[h],rhod,ff,grad[kk+i]);
      //printf("par grad: %f\n",grad[kk+i]);
    kk++;}
      h++;
    }
    
    Free(parC);
    return;
}


void Grad_Pair_Gauss_v(double rho, int *flag, double *gradcor, double *grad,
                     int *npar, double *par, double u, double v,
                     int *cormod,double lags,double lagt,int *flagcor, int *nparcT,double *parcor)
{
    // Initialization variables:
    int j=0;
    double mean=par[0],nugget=par[1],sill=par[2], *parC;
    
    parC=(double *) Calloc(nparcT[0],double);
    //printf("------- npar[0]: %d nparcT[0]: %d\n",npar[0],nparcT[0]); 4, 6
    for(j=0;j<nparcT[0];j++) parC[j]=parcor[j];
   
    double a=sill,b=sill*rho*(1-nugget),pa=a*a,pb=b*b;
    double c=-pa+pb,d=pa+pb,k=1/(c*c);
    double C=0.0,L=0.0,R=0.0, delta=0.0,rhod,zi,zj;
    double pu=0.0, pv=0.0;
    int h=0, i=0,  kk;
    //defines useful quantities:
    zi=u-mean;
    zj=v-mean;
    pu=pow(zi,2);
    pv=pow(zj,2);
    R=pu+pv;L=zi*zj;
    
    double ff=log_biv_Norm(rho*(1-nugget),u,v, mean, mean,sill,nugget);
    // Derivatives  respect with the mean
    if(flag[0]==1)
    {
        grad[i]= (((zi*a - zj*b)) + ((zj*a - zi*b)))/(c);i++;
        
        
    }
    // Derivative  respect with the nugget
    if(flag[1]==1)
    {
    delta=SQE1*nugget;
    grad[i]=(log_biv_Norm(rho*(1-(nugget+delta)),u,v, mean, mean,sill,0) -ff)/(delta);
        i++;
        //printf("nugget grad: %f\n",grad[i]);
       // grad[i]=0.5*k*(R*d-L*4*b*a-2*a*(pa-pb));i++;
    }
    // Derivative respect with the sill
    if(flag[2]==1)
    {
        delta=SQE1*sill;
        grad[i]=(log_biv_Norm(rho*(1-(nugget)),u,v, mean, mean,sill+delta,0) -ff)/(delta);
        //printf("sill grad: %f log: %f ff: %f SQE: %f nugget %f delta: %f\n",grad[i],log_biv_Norm(rho,u,v, mean, mean,sill+delta,nugget),ff,SQE, nugget,delta);
        //printf("sill grad: %f\n",grad[i]);
        i++;
    }
    // Derivatives with respect to the correlation parameters
    
    /*
    h=0;
    C=-k*sill*(R*a*b-L*d+b*c);
    for(j=i;j<*npar;j++){
        //printf("cr grad: %f \n",grad[i]);
        grad[j]=C*gradcor[h];h++;
        
    }*/
     
    h=0;
     kk=0;
  for(j=i;j<(i+*nparcT);j++) {
     
  if(flagcor[h]==1){
       delta=SQE1*parcor[h];
       parC[h]=parcor[h]+delta;
       rhod=CorFct(cormod,lags,lagt,parC);
      
      //printf("kk+i: %d h: %d rhod: %f lags: %f lagt: %f parC[h] %f parC[0] %f parC[1] %f parC[2] %f parC[3] %f parC[4] %f parC[5] %f parC[6] %f \n",kk+i,h,rhod,lags,lagt,parC[h],parC[0],parC[1],parC[2],parC[3],parC[4],parC[5],parC[6]);
      
     /* if(ISNAN(rhod))
      {
          printf("h: %d rhod: %f lags: %f lagt: %f parC[h] %f parC[0] %f parC[1] %f parC[2] %f parC[3] %f parC[4] %f parC[5] %f parC[6] %f \n",h,rhod,lags,lagt,parC[h],parC[0],parC[1],parC[2],parC[3],parC[4],parC[5],parC[6]);
      }*/
       
    grad[kk+i]=(log_biv_Norm(rhod*(1-nugget),u,v, mean, mean,sill,0) - ff)/delta;
      //printf("h: %d delta: %f parC[h]: %f rhod %f ff %f grad[kk+i]: %f\n",h, delta,parC[h],rhod,ff,grad[kk+i]);
      //printf("par grad: %f\n",grad[kk+i]);
    kk++;}
      h++;
    }
    
    Free(parC);
    return;
}


/********************************************************************/
/************ derivative of correlation function ****************/
/********************************************************************/

// Derivatives with respect to the separable parameter of the Gneiting correlation model:
double DGneiting_sep(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
    double a=0,arg=0,rho=0;
    arg=1+pow(u/scale_t, power_t);
    rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
    if(arg) a=0.5*pow(h/scale_s, power_s)*pow(arg,-0.5*sep*power_s)*power_s*rho*log(arg);
    return(a);
}
// Derivatives with respect to spatial power of the Gneiting correlation model:
double DGneiting_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
    double arg=0.0,rho=0.0,a=0.0;
    arg=1+pow(u/scale_t, power_t);
    rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
    if(h && arg)
    {
        a=(-pow(h/scale_s, power_s)*pow(arg,-0.5*sep*power_s)*log(h/scale_s) + 0.5*pow(h/scale_s, power_s)*pow(arg,-0.5*sep*power_s)*sep*log(arg) )*rho;
    }
    return(a);
}

// Derivatives with respect to temporal power of the Gneiting correlation model:
double DGneiting_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep)
{
    double arg=0.0,rho=0.0,a=0.0;
    arg=1+pow(u/scale_t, power_t);
    rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
    if(u)
    {
        a=( -pow(u/scale_t, power_t)*log(u/scale_t)*rho+0.5*rho*pow(h/scale_s, power_s)*power_s*sep*log(u/scale_t)*pow(arg,-0.5*sep*power_s))/arg;
    }
    return(a) ;
}

// Derivatives with respect to the temporal scale parameter of the Gneiting correlation model:
double DGneiting_sc_t(double h,double u, double power_s,double power_t,
                      double scale_s,double scale_t,double sep)
{
    double arg=0.0,rho=0.0,a=0.0;
    arg=1+pow(u/scale_t, power_t);
    rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
    if (arg) a= (power_t*pow(u/scale_t, power_t)*rho)/(scale_t*pow(arg,2))-( 0.5*power_t*power_s*sep*rho*pow(h/scale_s, power_s)*pow(u/scale_t, power_t)*pow(arg,-0.5*sep*power_s-2))/scale_t;
    return(a);
}
// Derivatives with respect to the spatial scale parameter of the Gneiting correlation model:
double DGneiting_sc_s(double h,double u,double power_s,double power_t,
                      double scale_s,double scale_t,double sep)
{
    double arg=0.0,rho=0.0,a=0.0;
    arg=1+pow(u/scale_t, power_t);
    rho=exp(-pow(h/scale_s, power_s)/(pow(arg, 0.5*sep*power_s)))/arg;
    a=(pow(h/scale_s, power_s)*power_s*rho*pow(arg,-0.5*power_s*sep-1))/scale_s;
    
    return (a);
}

// Derivatives with respect to scale of the Stable correlation model:
double DStabSc(double lag, double power, double scale, double rho)
{
    if(lag) return rho*power*pow(lag/scale,power)/scale;
    else return 0.0;
}
/********************************************************************/
/************** correlation  models **********************************/
/********************************************************************/
double CorFunBohman(double lag,double scale)
{
    double rho=0.0,x=0;
    x=lag/scale;
    if(x<=1)
    {
        if (x>0) rho=(1-x)*(sin(2*M_PI*x)/(2*M_PI*x))+(1-cos(2*M_PI*x))/(2*M_PI*M_PI*x);
        else   rho=1;
    }
    else rho=0;
    return rho;
}

// Stable class of correlation models:
double CorFunStable(double lag, double power, double scale)
{
    double rho=0.0;
    // Computes the correlation:
    rho=exp(-pow(lag/scale,power));
    return rho;
}


double CorFct(int *cormod, double h, double u, double *par)
{
   
    double arg=0.0,  power_s=0.0, power_t=0.0;
    double rho=0.0, sep=0,  scale_s=0.0, scale_t=0;
    
    switch(*cormod) // Correlation functions are in alphabetical order
    {
            // START non-separable correlation functions:
        case 1:// Double exp:
            scale_s=par[0];
            scale_t=par[1];
            rho=CorFunStable(h,1,scale_s)*CorFunStable(u,1,scale_t);
            break;
        case 2:   //Gneiting correlation model as in (14) Gneitint (2002) with tau=1
            power_s=par[0];
            power_t=par[1];
            scale_s=par[2];
            scale_t=par[3];
            sep=par[4];
            arg=1+pow(u/scale_t, power_t);
            rho=exp(-(pow(h/scale_s, power_s))*pow(arg, -0.5*sep*power_s))/pow(arg,1);
            break;
            
        case 3:   //WendLand Time correlation model
            rho=wen_time(par,h,u);
            break;
    }
    //printf("rho: %f \n",rho);
    return rho;
}


double CorFct_v(int *cormod, double h, double u, double *par)
{
   
    double arg=0.0,  power_s=0.0, power_t=0.0;
    double rho=0.0, sep=0,  scale_s=0.0, scale_t=0;
    
    switch(*cormod) // Correlation functions are in alphabetical order
    {
            // START non-separable correlation functions:
        case 1:// Double exp:
            scale_s=par[0];
            scale_t=par[1];
            rho=CorFunStable(h,1,scale_s)*CorFunStable(u,1,scale_t);
            break;
        case 2:   //Gneiting correlation model as in (14) Gneitint (2002) with tau=1
            power_s=par[0];
            power_t=par[1];
            scale_s=par[2];
            scale_t=par[3];
            sep=par[4];
            arg=1+pow(u/scale_t, power_t);
            rho=exp(-(pow(h/scale_s, power_s))*pow(arg, -0.5*sep*power_s))/pow(arg,1);
            break;
            
        case 3:   //WendLand Time correlation model
            rho=wen_time(par,h,u);
            break;
    }
    //printf("rho: %f \n",rho);
    return rho;
}


// cdf of  a bivariate Gausssian distribution
void CorFct_call(int *cormod, double *h, double *u, double *par,double *res)
{
    
    *res=   CorFct(cormod, h[0], u[0], par);
    //Rprintf("rho: %f\n",*res);
    //*res=   cdf_norm(0,0,.5,1);
}
/********************************************************************/
/********************************************************************/
/********************************************************************/

//  Derivatives with respect ot the correlations parameters:
void GradCorrFct_v(double rho, int *cormod, int *flag,
                 double *grad, double h, double u,double *par)
{
    int i=0;
    double power_s=0.0, power_t=0.0;
    double scale_s=0.0, scale_t=0.0, sep=0.0;
    
    switch(*cormod)// Correlation functions are in alphabetical order
    {
            //spatial gradients of correlations:
        case 1://Double Exponentil
            scale_s=par[0];
            scale_t=par[1];
            if(flag[0]==1){//spatial-scale parameter
                grad[i]=DStabSc(h,1,scale_s,rho);i++;}
            //temporal-scale parameter
            if(flag[1]==1) grad[i]=DStabSc(u,1,scale_t,rho);
            break;
        case 2://Gneiting spatio-temporal correlation
            power_s=par[0];
            power_t=par[1];
            scale_s=par[2];
            scale_t=par[3];
            sep=par[4];
            
            if(flag[0]==1){//spatial-power parameter
                grad[i]=DGneiting_pw_s(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
            if(flag[1]==1){//temporal-power parameter
                grad[i]=DGneiting_pw_t(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
            if(flag[2]==1){//spatial-scale parameter
                grad[i]=DGneiting_sc_s(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
            if(flag[3]==1){//temporal-scale parameter
                grad[i]=DGneiting_sc_t(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
            //separable parameter
            if(flag[4]==1)
                grad[i]=DGneiting_sep(h,u,power_s,power_t,scale_s,scale_t,sep);
            break;
        case 3://WendLand time

            //0:power2_s 1:power_s=2  2:power2_t   3:scale_s     4:scale_t   5:sep 6:smooth
            if(flag[0]==1){//power2_s parameter
                grad[i]=deri_R_power_wen_time(par, h,u);i++;
            }
            if(flag[1]==1){//power_s parameter
                grad[i]=0;i++;}
            if(flag[2]==1){//power2_t parameter
                grad[i]=deri_R_power_t_wen_time_v(par, h,u);i++;}
            if(flag[3]==1){//scale_s parameter
                //grad[i] = 0.1;
                grad[i]=deri_scale_s_wen_time_v(par, h,u);i++;}
            if(flag[4]==1){//scale_t parameter
                //grad[i] = 0.1;
                grad[i]=deri_scale_t_wen_time_v(par, h,u);i++;}
            if(flag[5]==1){//sep parameter
                //grad[i] = 0.1;i++;}
                grad[i]=deri_sep_wen_time_v(par, h,u);i++;}
            //smooth parameter
            if(flag[6]==1)
            {
                //grad[i] = 0.1;
                grad[i]=deri_smooth_wen_time_v(par, h,u);
            }
            break;
    }
    
}


//  Derivatives with respect ot the correlations parameters:
void GradCorrFct(double rho, int *cormod, int *flag,
                 double *grad, double h, double u,double *par)
{
    int i=0;
    double power_s=0.0, power_t=0.0;
    double scale_s=0.0, scale_t=0.0, sep=0.0;
    
    switch(*cormod)// Correlation functions are in alphabetical order
    {
            //spatial gradients of correlations:
        case 1://Double Exponentil
            scale_s=par[0];
            scale_t=par[1];
            if(flag[0]==1){//spatial-scale parameter
                grad[i]=DStabSc(h,1,scale_s,rho);i++;}
            //temporal-scale parameter
            if(flag[1]==1) grad[i]=DStabSc(u,1,scale_t,rho);
            break;
        case 2://Gneiting spatio-temporal correlation
            power_s=par[0];
            power_t=par[1];
            scale_s=par[2];
            scale_t=par[3];
            sep=par[4];
            
            if(flag[0]==1){//spatial-power parameter
                grad[i]=DGneiting_pw_s(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
            if(flag[1]==1){//temporal-power parameter
                grad[i]=DGneiting_pw_t(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
            if(flag[2]==1){//spatial-scale parameter
                grad[i]=DGneiting_sc_s(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
            if(flag[3]==1){//temporal-scale parameter
                grad[i]=DGneiting_sc_t(h,u,power_s,power_t,scale_s,scale_t,sep);i++;}
            //separable parameter
            if(flag[4]==1)
                grad[i]=DGneiting_sep(h,u,power_s,power_t,scale_s,scale_t,sep);
            break;
        case 3://WendLand time

            //0:power2_s 1:power_s=2  2:power2_t   3:scale_s     4:scale_t   5:sep 6:smooth
            if(flag[0]==1){//power2_s parameter
                grad[i]=deri_R_power_wen_time(par, h,u);i++;
            }
            if(flag[1]==1){//power_s parameter
                grad[i]=0;i++;}
            if(flag[2]==1){//power2_t parameter
                grad[i]=deri_R_power_t_wen_time(par, h,u);i++;}
            if(flag[3]==1){//scale_s parameter
                //grad[i] = 0.1;
                grad[i]=deri_scale_s_wen_time(par, h,u);i++;}
            if(flag[4]==1){//scale_t parameter
                //grad[i] = 0.1;
                grad[i]=deri_scale_t_wen_time(par, h,u);i++;}
            if(flag[5]==1){//sep parameter
                //grad[i] = 0.1;i++;}
                grad[i]=deri_sep_wen_time(par, h,u);i++;}
            //smooth parameter
            if(flag[6]==1)
            {
                //grad[i] = 0.1;
                grad[i]=deri_smooth_wen_time(par, h,u);
            }
            break;
    }
    
}

// Computes the Geodesic distance between to coordinates:
double Dist_geodesic(double loni, double lati, double lonj, double latj)
{
    double ai, bi, aj, bj, val;
    val = 0.0;
    if (loni == lonj && lati == latj) return val;
    ai = (lati)*M_PI/180;
    bi = (loni)*M_PI/180;
    aj = (latj)*M_PI/180;
    bj = (lonj)*M_PI/180;
    val = sin(ai) * sin(aj) + cos(ai) * cos(aj) * cos(bi - bj);
    val = acos(val) *  REARTH;
    return val;
}


/**************** return a bivariate vector with maximum and minimum of a vector *********/
void Range(double *x, double *ran, int size)
{
    int i=0;
    
    ran[0] = x[0];
    ran[1] = x[0];
    
    for(i = 1; i < size; i++)
    {
        ran[0] = fmin(ran[0], x[i]);
        ran[1] = fmax(ran[1], x[i]);
    }
    return;
}
/**********************************************************************************/
int is_equal(double val1, double val2)
{
    return fabs(val1-val2)<MAXERR;
}


void SeqStep(double *x, int len, double step, double *res)
{
    int i=0;
    res[0]=x[0];
    for(i=0;i<len;i++) {
        res[i+1]=res[i]+step;
    
    }
    return;
}
/****************************************************************************************************************/
/****************************************************************************************************************/
/****************************************************************************************************************/
// subsampling in space
void SetSampling_s(int ncoord,int ntime,double *coordx, double *coordy, double *data, int *npts,double *scoordx, double *scoordy, double *sdata, double xmax,double xmin, double ymax,double ymin)
{
    int i=0, j=0, f=0,h=0;
    for(i=0;i<ncoord;i++)
    {
        if((xmin<coordx[i]||is_equal(xmin,coordx[i]))&&
           (xmax>coordx[i]||is_equal(xmax,coordx[i]))&&
           (ymin<coordy[i]||is_equal(ymin,coordy[i]))&&
           (ymax>coordy[i]||is_equal(ymax,coordy[i])))
        {
            scoordx[j]=coordx[i];
            scoordy[j]=coordy[i];
           
            for(h=0;h<ntime;h++)
            {
                sdata[f]=data[h+ntime*i];f++;
            }
            j++;
        }
    }
    *npts = j;
    return;
}
// subsampling in time
void SetSampling_t(double *data,double *sdata,int ncoord,int ntime,int wint,int k)
{
    int i=0,j=0,p=0;
    for(i=0;i<(ncoord);i++)
        for(j=(k+(ntime*i));j<(k+wint+(ntime*i));j++)
        {
            sdata[p]=data[j];p++;
        }
    
    return;
}
