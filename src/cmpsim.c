#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
/*#include <omp.h>*/

#define NITER 10000

double logsumexp1(double a, double b)
{
  double r;
  if(a > b) r = a + log(1 + exp(b-a));
  else r = b + log(1 + exp(a-b));
  return(r);
}

/* log Z function*/

double compoisson_z_ln1(double lamb, double nusb)
{
  double tol = 1e-8;
  double lambda_ln;
  double z = log(1+lamb), a, r;
  int i=2;
  if(lamb <= 0) {
    nusb = 0;
  }
  if(nusb == 0) {
    /* special case */
    /* needed when *lambda is close to 1 */
    z = -log(1-lamb);
  }
  a = lambda_ln = log(lamb);
  for(;i<NITER;i++) {
    double e = lambda_ln - nusb*log(i);
    a = a+e;
    z = logsumexp1(z,a);
    if(e < 0) {
      /* truncation error */
      r = exp(a-z)/(1-exp(e));
      /* truncation error in log domain */
      /* log(z+r)-log(z) = log(1 + r/z) =approx r/z */
      if(r < tol) break;
    }
  }
  return z;
}


/* Simulation*/

double cmpsim(double lam,double nus)
{
  double u;
  double cnt,mf;
  /*time_t t;*/
  /*  */
    cnt=0;
/* Intializes random number generator */
   /*srand((unsigned) time(&t));*/
    u=runif(0,1); /*/(double)RAND_MAX;*/
/*printf("l: %1f \n",lam);
printf("n: %1f \n",nus);*/
    double z_ln=compoisson_z_ln1(lam,nus);
/*printf("z: %1f \n",z_ln);*/
    double proby=1/exp(z_ln);
double pres=1/exp(z_ln);
    double lam_ln=log(lam);
/*printf("prob: %1f \n",proby);*/
while(proby<u){
      cnt=cnt+1;
/*printf("%1f \n",cnt);*/
mf=lam_ln-nus*log((double)cnt);
if(pres==0)
{
 pres=exp(mf-z_ln);
 proby=proby+pres;
}else
{
pres=pres*exp(mf);
proby=proby+pres;
/*printf("probyi:%1f \n",proby);*/
 }
}
  return cnt;
}


void  cmpsim_all(double *lambda,double *nu,int *n,double *y)
{
int i=0;
int n1=*n;
GetRNGstate();
/*int iCPU=omp_get_num_procs();
omp_set_num_threads(iCPU);
#pragma omp parallel for */
for(i=0;i<n1;i++){
y[i]=cmpsim(lambda[i],nu[i]);
}

PutRNGstate();

}



