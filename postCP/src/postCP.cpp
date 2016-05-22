#define INF 1e300
//#include "postCP.h"

#include <R.h>
#include <Rmath.h>
#include <stddef.h>
#include <vector>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

#include <ctime>

#define INF 1e300


extern "C"{
using std::vector;
using std::string;
using std::cout;
using std::setprecision;
using std::endl;

// #include <R.h>

using namespace std;

//vector<double> data;
vector<int> seg;
vector<double> mu;


int n;
int J;
int nseg;
int mout;

int normal;
double sigma;
double ci;
int probs;

int priortype;
int nsamples;
int verbose,debug;
double priortemp;
int levelbased;

double diffclock(clock_t clock1,clock_t clock2) {
	double diffticks=clock1-clock2;
	double diffs=(diffticks)/CLOCKS_PER_SEC;
	return diffs;
}

// returns log(exp(la)+exp(lb))
double lsum(double la, double lb) {
  if (la>lb){
    if (lb>(-INF)) return la+log1p(exp(lb-la));
    else return la;
    }
  else{
    if (la>(-INF)) return lb+log1p(exp(la-lb));    
    else return lb;
  }
}

void forward(double* data, double* lforward, double* lprob, double* rprior) {
  vector<double> ldist(J);

  double aux1=-0.5*log(2*M_PI)-log(sigma);
  double aux2=-0.5;
  //double aux3=log(0.5);
 
 // bool lower;
if (debug) Rprintf("normal?f %i\n",normal);
 for (int j=0; j<J; j++) for (int i=0; i<n; i++) lforward[n*j+i]=-INF;

  if (!probs){
  if (normal) lforward[0]=aux1+aux2*(data[0]-mu[0])*(data[0]-mu[0])/sigma/sigma;
  else lforward[0]=data[0]*log(mu[0])-mu[0]-lgamma(data[0]+1);} else{
     lforward[0]=lprob[0];
  }
//  for (int j=1; j<J; j++) lforward[n*j]=-INF;
  if (verbose) Rprintf("\nCalculating forward matrix...\n");     
//  int a=0;
  for (int i=1; i<n; i++) {
    if (verbose) {
      if ((i+1)%10000==0)   Rprintf(".");
    }
    // precompute ldist for data[i]
    for (int j=0; j<J; j++) {
      if (!probs){
        if (normal) ldist[j]=aux1+aux2*(data[i]-mu[j])*(data[i]-mu[j])/sigma/sigma;
        else ldist[j]=data[i]*log(mu[j])-mu[j]-lgamma(data[i]+1);
      } else ldist[j]=lprob[n*j+i];
         
     }  
    // recursion step
    if (priortype==2) priortemp=rprior[i-1];
    if (priortype==3) priortemp=rprior[0];

    lforward[i]=lforward[i-1]+log(1-priortemp)+ldist[0];

    for (int j=1; j<J; j++) {
    if (priortype==2) priortemp=rprior[i-1];
    if (priortype==3) priortemp=rprior[j];

      lforward[n*j+i]=lsum(lforward[n*(j-1)+i-1]+log(priortemp),lforward[(n*j)+i-1]+log(1-priortemp))+ldist[j];
    }
  }  
}



void backward(double* data, double* lbackward, double* lforward, double* lprob, double* cp, double* rprior) {
  vector<double> ldist(J);

  double aux1=-0.5*log(2*M_PI)-log(sigma);
  double aux2=-0.5;
  double priortemp1;
 // double aux3=log(0.5);

 // bool upper;    
  if (debug) Rprintf("normal?b %i\n",normal);
 for (int j=0; j<J; j++) for (int i=0; i<n; i++) lbackward[n*j+i]=-INF;

//  int a=J-2;
  lbackward[n*J-1]=0.0;
//  for (int j=0; j<(J-1); j++) lbackward[n*j+n-1]=-INF;
  if (verbose) Rprintf("\nCalculating backward matrix...\n");      
  for (int i=n-2; i>=0; i--) {
   if (verbose) {
      if ((i+1)%10000==0)  Rprintf("."); 
    }
    // precompute ldist for data[i+1]
    for (int j=0; j<J; j++) {
      if (!probs){
        if (normal) ldist[j]=aux1+aux2*(data[i+1]-mu[j])*(data[i+1]-mu[j])/sigma/sigma;
        else ldist[j]=data[i+1]*log(mu[j])-mu[j]-lgamma(data[i+1]+1);} else{
           ldist[j]=lprob[n*j+i+1];
       }
    }
    // recursion step
    if (priortype==2) priortemp=rprior[i];
    if (priortype==3) priortemp=rprior[J-1];

     lbackward[n*(J-1)+i]=lbackward[n*(J-1)+i+1]+log(1-priortemp)+ldist[J-1];
    for (int j=J-2; j>=0; j--) {
    if (priortype==2) priortemp=rprior[i];
    if (priortype==3) {
      priortemp=rprior[j];
      priortemp1=rprior[j+1];
    } else priortemp1=priortemp;
 
    lbackward[n*j+i]=lsum(lbackward[n*j+i+1]+ldist[j]+log(1-priortemp),lbackward[n*(j+1)+i+1]+ldist[j+1]+log(priortemp1)); 
    cp[(n-1)*j+i]=exp(lforward[n*j+i]+lbackward[n*(j+1)+i+1]+ldist[j+1]+log(priortemp1)-lforward[n*J-1]-lbackward[n*J-1]);
    if (isnan(cp[(n-1)*j+i])) cp[(n-1)*j+i]=0;
    }
  }  

}



void forward_lev(double* data, double* lforward, double* lprob, double* rprior) {
  vector<double> ldist(J);

  double aux1=-0.5*log(2*M_PI)-log(sigma);
  double aux2=-0.5;
  //double aux3=log(0.5);
 
 // bool lower;
if (debug) Rprintf("normal? %i\n",normal);
 for (int j=0; j<J; j++) for (int k=0; k<nseg; k++) for (int i=0; i<n; i++) lforward[n*nseg*j+n*k+i]=-INF;

  if (!probs){
    if (normal) lforward[0]=aux1+aux2*(data[0]-mu[0])*(data[0]-mu[0])/sigma/sigma;
    else lforward[0]=data[0]*log(mu[0])-mu[0]-lgamma(data[0]+1);
  } else{
    lforward[0]=lprob[0];
  }

  if (verbose) Rprintf("\nCalculating forward matrix...\n");     
//  int a=0;
  for (int i=1; i<n; i++) {
    if (verbose) {
      if ((i+1)%10000==0)   Rprintf(".");
    }
    // precompute ldist for data[i]
    for (int j=0; j<J; j++) {
      if (!probs){
        if (normal) ldist[j]=aux1+aux2*(data[i]-mu[j])*(data[i]-mu[j])/sigma/sigma;
        else ldist[j]=data[i]*log(mu[j])-mu[j]-lgamma(data[i]+1);} else{
       ldist[j]=lprob[n*j+i];
      }
    }
    // recursion step
    if (priortype==2) priortemp=rprior[i-1];
    if (priortype==3) priortemp=rprior[0];
    lforward[i]=lforward[i-1]+log(1-priortemp)+ldist[0];
    for (int j=1; j<J; j++) lforward[n*nseg*j+i]=-INF;      

    for (int k=1; k<nseg; k++) {
     // lforward[n*k+i]=lforward[n*k+i-1]+log(1-priortemp)+ldist[j];
      for (int j=0; j<J; j++) {
        if (priortype==2) priortemp=rprior[i-1];
        if (priortype==3) priortemp=rprior[k];    
        lforward[n*nseg*j+n*k+i]=lforward[n*nseg*j+n*k+i-1]+log(1-priortemp);      
          for (int jj=0; jj<J; jj++) {
             if (jj != j) lforward[n*nseg*j+n*k+i]=lsum(lforward[n*nseg*j+n*k+i],lforward[n*nseg*jj+n*(k-1)+i-1]+log(priortemp)-log(J-1.0));
         }
        lforward[n*nseg*j+n*k+i]=lforward[n*nseg*j+n*k+i]+ldist[j];
      }
    }
  }  
}

void backward_lev(double* data, double* lbackward, double* lforward, double* lprob, double* cp, double* rprior) {
  vector<double> ldist(J);

  double aux1=-0.5*log(2*M_PI)-log(sigma);
  double aux2=-0.5;
  double priortemp1;
 // double aux3=log(0.5);

 // bool upper;    
  if (debug) Rprintf("normal? %i\n",normal);
  for (int j=0; j<J; j++) for (int k=0; k<nseg; k++) for (int i=0; i<n; i++) lbackward[n*nseg*j+n*k+i]=-INF;
  for (int j=0; j<J; j++) lbackward[n*nseg*j+n*(nseg-1)+n-1]=0;

  double evidence1=-INF;  
  for (int j=0; j<J; j++) evidence1=lsum(evidence1,lforward[n*nseg*j+n*(nseg-1)+n-1]+lbackward[n*nseg*j+n*(nseg-1)+n-1]);
  if (verbose) Rprintf("\nCalculating backward matrix...\n");      
  for (int i=n-2; i>=0; i--) {
    if (verbose) {
      if ((i+1)%10000==0)  Rprintf("."); 
    }
    // precompute ldist for data[i+1]
    for (int j=0; j<J; j++) {
      if (!probs){
        if (normal) ldist[j]=aux1+aux2*(data[i+1]-mu[j])*(data[i+1]-mu[j])/sigma/sigma;
        else ldist[j]=data[i+1]*log(mu[j])-mu[j]-lgamma(data[i+1]+1);} else{
         ldist[j]=lprob[n*j+i+1];
      }
    }
    // recursion step
    if (priortype==2) priortemp=rprior[i];
    if (priortype==3) priortemp=rprior[nseg-1];
    for (int j=J-1; j>=0; j--) lbackward[n*nseg*j+n*(nseg-1)+i]=lbackward[n*nseg*j+n*(nseg-1)+i+1]+ldist[j]+log(1-priortemp);   

//   Rprintf("%0.8f%0.8f.%0.8f.\n",data[i+1],mu[j],ldist[j]);}
    for (int k=(nseg-2); k>=0; k--) {
      cp[(n-1)*k+i]=0;

      for (int j=J-1; j>=0; j--) {
        if (priortype==2) priortemp=rprior[i];
        if (priortype==3) {
         priortemp=rprior[k];
         priortemp1=rprior[k+1];
      } else priortemp1=priortemp;
        lbackward[n*nseg*j+n*k+i]=lbackward[n*nseg*j+n*k+i+1]+ldist[j]+log(1-priortemp);
        for (int jj=J-1; jj>=0; jj--) {
          if (jj != j) {
            lbackward[n*nseg*j+n*k+i]=lsum(lbackward[n*nseg*j+n*k+i],lbackward[n*nseg*jj+n*(k+1)+i+1]+ldist[jj]+log(priortemp1)-log(J-1.0)); 
            cp[(n-1)*k+i]=cp[(n-1)*k+i]+exp(lforward[n*nseg*j+n*k+i]+lbackward[n*nseg*jj+n*(k+1)+i+1]+ldist[jj]+log(priortemp1)-log(J-1.0)-evidence1);
          }
        }
      }
    }  
  }
  double evidence2=-INF;
 for (int j=0; j<J; j++) evidence2=lsum(evidence2,lforward[n*nseg*j]+lbackward[n*nseg*j]);
   if (verbose) Rprintf("\n log(pevidence)=%.8f=%.8f\n",evidence2,evidence1);

}



void postCPsample(double* data, double* rmu, double *rsigma, double* lprob, int *rprobs,  double* rlbackward, double* cpvector, double* levvector, int *rnsamples, int *nn, int *JJ, int *rnseg, int *rmodel, double* rprior, int *rpriortype, int *rlevelbased, int *rverbose, int *rdebug) {


int nsamples=*rnsamples;
int n=*nn;
int J=*JJ;
//int jj;
if (debug) Rprintf("start");
int nseg=*rnseg;
if (debug) Rprintf(" sampling");

int verbose=*rverbose;
int debug=*rdebug;
int priortype=*rpriortype;

if (debug){ 
Rprintf("\nGenerating %i change-points,",nseg-1);
Rprintf(" %i replicates\n",nsamples);
}
probs=*rprobs;
normal=0;
levelbased=*rlevelbased;
if (*rmodel==2) normal=1;
double sigma=*rsigma;
 GetRNGstate(); 

    if (priortype==1) priortemp=rprior[0];

for (int r=0; r<nsamples; r++) {
  // int cp=0;
   int i=0;
  
   double rand;
   double aux=0;
  
   double aux1=-0.5*log(2*M_PI)-log(sigma);
   double aux2=-0.5;
//   double aux3=log(0.5);
  
   if (verbose&((r+1)%500==0)) Rprintf("\nGenerating replicate %i\n",r+1);
   int j=0;
   int jj=1;
   rand=runif(0,1.0);
   for (int k=0;k<nseg-1;k++){
     aux=0;
     while (true) {
   //  if (debug)  Rprintf("i:%i,k:%i,n:%i,nseg:%i",i,k,n,nseg);   
     if (priortype==2) priortemp=rprior[i];
     if (priortype==3) priortemp=rprior[k];
       if (!levelbased){
         if (!probs){
             if (debug)  Rprintf("%i,%i:%.8f\n",i,k,aux);
          if (normal) aux=rlbackward[n*(k+1)+i+1]-rlbackward[n*k+i]+log(priortemp)+aux1+aux2*(data[i+1]-rmu[k+1])*(data[i+1]-rmu[k+1])/sigma;
          else aux=rlbackward[n*(k+1)+i+1]-rlbackward[n*k+i]+log(priortemp)+data[i+1]*log(rmu[k+1])-rmu[k+1]-lgamma(data[i+1]+1);
         } else aux=rlbackward[n*(k+1)+i+1]-rlbackward[n*k+i]+log(priortemp)+lprob[n*(k+1)+i+1]; 
         rand=log(runif(0,1.0));
         i++;
         if (rand<aux)
         break;} // if !levelbased
       else{
         if (jj==0) rand=runif(0,1.0);
         if (jj!=j){ 
          if (!probs){      
  	  if (normal) aux=aux+exp(rlbackward[n*nseg*jj+n*(k+1)+i+1]-rlbackward[n*nseg*j+n*k+i]+log(priortemp)-log(J-1.0)+aux1+aux2*(data[i+1]-rmu[jj])*(data[i+1]-rmu[jj])/sigma);
          else aux=aux+exp(rlbackward[n*nseg*jj+n*(k+1)+i+1]-rlbackward[n*nseg*j+n*k+i]+log(priortemp)-log(J-1.0)+data[i+1]*log(rmu[jj])-rmu[jj]-lgamma(data[i+1]+1));
	  } else aux=aux+exp(rlbackward[n*nseg*jj+n*(k+1)+i+1]-rlbackward[n*nseg*j+n*k+i]+log(priortemp)-log(J-1.0)+lprob[n*jj+i+1]);
         } // jj!=j
   	 if (rand<aux){
           j=jj;
           jj=0;
           i++;          
           break;
         }  // rand<aux         
         if (jj==J-1) {i++;jj=0;aux=0;}
         else jj++;
         
   
       } // else levelbased            
     } // end while loop
     cpvector[nsamples*k+r]=i;     
     if (levelbased) {levvector[nsamples*k+r]=j+1;
	 if (debug)  Rprintf("cp:%i,l:%i",i,j+1);}
   } // end for loop
 } // end r loop
PutRNGstate();
//   return;
}




void viterbi(double* data, double* lprob, double* bestcp, double* rprior){
  if (debug) Rprintf("viterbi calculations:%i x %i!!,\n",n,J);
  if (debug) Rprintf("n%i J%i nseg%i",n,J,nseg); 
  double* lviterbi=NULL;
  lviterbi=new double[n*J];
  int* vitstate=NULL;
  vitstate=new int[n*J];
  double aux1=-0.5*log(2*M_PI)-log(sigma);
  double aux2=-0.5;
  double ldist;
  if (debug) Rprintf("%.8f\n",lprob[0]);
  // bool lower;
 int change; 
  if (verbose) Rprintf("normal?v %i\n",normal);
  for (int j=0;j<J;j++)   for (int i=0;i<n;i++) lviterbi[n*j+i]=-INF;

  if (!probs){
    if (normal) ldist=aux1+aux2*(data[0]-mu[0])*(data[0]-mu[0])/sigma/sigma;
    else ldist=data[0]*log(mu[0])-mu[0]-lgamma(data[0]+1);} else{
      ldist=lprob[0];
    }
    lviterbi[0]=ldist;
; 

  if (verbose) Rprintf("\nApplying Viterbi algorithm ...\n");      

         vitstate[0]=0; 
  for (int i=1; i<n; i++) {
     if (priortype==2) priortemp=rprior[i-1];
     if (priortype==3) priortemp=rprior[0];
     if (!probs){
       if (normal) ldist=aux1+aux2*(data[i]-mu[0])*(data[i]-mu[0])/sigma/sigma;
       else ldist=data[i]*log(mu[0])-mu[0]-lgamma(data[i]+1);} else{
            ldist=lprob[i];
      }
         lviterbi[i]=log(1-priortemp)+lviterbi[i-1]+ldist;
         vitstate[i]=0;     
   if (verbose)   if ((i+1)%10000==0)  Rprintf("."); 
     for (int j=1;j<J;j++){ 
       if (priortype==2) priortemp=rprior[i-1];
       if (priortype==3) priortemp=rprior[j];
       if (!probs){
         if (normal) ldist=aux1+aux2*(data[i]-mu[j])*(data[i]-mu[j])/sigma/sigma;
         else ldist=data[i]*log(mu[j])-mu[j]-lgamma(data[i]+1);}else{
            ldist=lprob[n*j+i];
         }
      // if (debug)Rprintf("%0.7f",ldist);
       change=(lviterbi[n*(j-1)+i-1]+log(priortemp))>(lviterbi[n*j+i-1]+log(1-priortemp));
       lviterbi[n*j+i]=ldist+(lviterbi[n*(j-1)+i-1]+log(priortemp))*change+(lviterbi[n*j+i-1]+log(1-priortemp))*(1-change);

       if (change){
         vitstate[n*j+i]=j-1; 
       } else{
         vitstate[n*j+i]=j; 
       }      
    if (debug) Rprintf("v:%i\n",vitstate[n*j+i]); 
     }
  } 
     int ii=n-1;
     int jj=J-1;
    delete [] lviterbi;
     while ((ii>=1)&(jj>=1)) {
       if ((vitstate[n*jj+ii]==jj-1)|(ii==1)){
          bestcp[jj-1]=ii;
          jj--;    
if (debug) Rprintf("cp:%i\n",ii);       
       }
     
       ii--;  
       if (ii==1) jj=1;
     }
   delete [] vitstate;
}


void viterbi_lev(double* data, double* lprob, double* bestcp, double* bestlevel, double* rprior){
  if (debug) Rprintf("viterbi calculations:%i x %i!!,\n",n,J);

  double* lviterbi=NULL;
  lviterbi=new double[n*nseg*J];
  int* vitstate=NULL;
  vitstate=new int[n*nseg*J];
  int* vitlevel=NULL;
  vitlevel=new int[n*nseg*J];
  int lastlevel=0;

  double maxvit;
  double tempvit;

  double aux1=-0.5*log(2*M_PI)-log(sigma);
  double aux2=-0.5;
  vector<double> ldist(J);
  if (debug) Rprintf("%.8f\n",lprob[0]);
  if (debug) Rprintf("THIS IS J!!!!!!!!!!!!!%i??",J);
  if (verbose) Rprintf("normal?v %i\n",normal);
  for (int j=1; j<J; j++) for (int k=1; k<nseg; k++) for (int i=0; i<n; i++)   lviterbi[n*nseg*j+n*k+i]=-INF;
  for (int k=1; k<nseg; k++) for (int i=0; i<n; i++)   lviterbi[n*k+i]=-INF; 
   for (int j=1; j<J; j++)  for (int i=0; i<n; i++)   lviterbi[n*nseg*j+i]=-INF;

  if (!probs){
    if (normal) ldist[0]=aux1+aux2*(data[0]-mu[0])*(data[0]-mu[0])/sigma/sigma;
    else ldist[0]=data[0]*log(mu[0])-mu[0]-lgamma(data[0]+1);
    } else{
      ldist[0]=lprob[0];
    }
  lviterbi[0]=ldist[0];


  if (verbose) Rprintf("\nApplying Viterbi algorithm ...\n");      
  //int  j=0;
  for (int k=0; k<nseg; k++) {vitstate[n*k]=0;
vitlevel[n*k]=0;}
//  int k=0;

  for (int i=1; i<n; i++) {
     if (priortype==2) priortemp=rprior[i-1];
     if (priortype==3) priortemp=rprior[0];
     if (!probs){
       if (normal) ldist[0]=aux1+aux2*(data[i]-mu[0])*(data[i]-mu[0])/sigma/sigma;
       else ldist[0]=data[i]*log(mu[0])-mu[0]-lgamma(data[i]+1);} else{
            ldist[0]=lprob[i];
      }
//   for (int j=0; j<J; j++)   
//    vitstate[i]=vitstate[i-1];     
    //  vitlevel[i]=vitlevel[i-1]; 
    if (verbose)   if ((i+1)%10000==0)  Rprintf("."); 
    for (int k=0;k<nseg;k++){ 
      // if (debug)Rprintf("%0.7f",ldist);
       if (priortype==2) priortemp=rprior[i-1];
       for (int j=0;j<J;j++){ 
         if (priortype==3) priortemp=rprior[0];
         maxvit=lviterbi[n*nseg*j+n*k+i-1]+log(1-priortemp);
  //       change=false;
         vitstate[n*nseg*j+n*k+i]=k; 
         vitlevel[n*nseg*j+n*k+i]=j; 
         if (k>0){
        
           if (priortype==3) priortemp=rprior[k-1];

  //         changevit-INF; 
           for (int jj=0;jj<J;jj++){  
             if (jj!=j){  
               tempvit=lviterbi[n*nseg*jj+n*(k-1)+i-1]+log(priortemp); // log(J-1) penalizes changes too much               
               if (tempvit>maxvit){
                 maxvit=tempvit;   
                 vitstate[n*nseg*j+n*k+i]=k-1; 
                 vitlevel[n*nseg*j+n*k+i]=jj; 
             //    change=true;
                 lastlevel=jj;
               }
             }             
           }
        } 

        if (!probs){
           if (normal) ldist[j]=aux1+aux2*(data[i]-mu[j])*(data[i]-mu[j])/sigma/sigma;
           else ldist[j]=data[i]*log(mu[j])-mu[j]-lgamma(data[i]+1);
           }else{
           ldist[j]=lprob[n*j+i];
         }
    
      lviterbi[n*nseg*j+n*k+i]=ldist[j]+maxvit;
      }
    }
  } 
     
    maxvit=-INF;
    for (int j=0;j<J;j++) {
        tempvit=lviterbi[n*nseg*j+n*(nseg-1)+n-1];
        if (debug) Rprintf("vj:%i,%0.7f,\n", j,tempvit);
        if (tempvit>maxvit) {
          maxvit=tempvit;
          lastlevel=j;
        }
     }   


     int ii=n-1;
     int kk=nseg-1;
     int    jj=lastlevel;
     bestlevel[nseg-1]=jj+1;
     delete [] lviterbi;
     while ((ii>=1)&(kk>=1)) {
       if ((vitstate[n*nseg*jj+n*kk+ii]==kk-1)|(ii==1)){
          jj=vitlevel[n*nseg*jj+n*kk+ii];
          bestcp[kk-1]=ii;
          bestlevel[kk-1]=jj+1; 
          kk--;          
       }
       ii--; 
     } 


if (debug){
    Rprintf("\n");
//    for (int j=0;j<J;j++){Rprintf("\n");  for (int i=0;i<n;i++) {Rprintf("\n"); for (int k=0;k<nseg;k++)  Rprintf("%.2f ",lviterbi[n*nseg*j+n*k+i]);}}

    Rprintf("\n");
    for (int i=0;i<n;i++) {Rprintf("\n"); for (int k=0;k<nseg;k++) Rprintf("%i ",vitstate[n*k+i]);}
    Rprintf("<-vitstates\n");
    for (int i=0;i<n;i++) {Rprintf("\n"); for (int j=0;j<J;j++) Rprintf("%i ",vitlevel[n*j+i]);}
    Rprintf("\n");
    Rprintf("<-vitlevels\n");
    for (int k=0;k<nseg-1;k++) {Rprintf("%.2f ",bestlevel[k]);}
}
   delete [] vitstate;
}


void viterbi2(double* data, double* lprob, int *rprobs,double* bestcp, double* bestlevel, double* rprior, int *rpriortype, int *nn, int *JJ, int *rnseg, double* rmu, double *rsigma,int *rmodel,int *level, int *rverbose, int *rdebug){
  verbose=*rverbose;
  debug=*rdebug;
   n=*nn; 
   J=*JJ;
   nseg=*rnseg;
   probs=*rprobs;
   int level2=*level;
   mu.resize(J);
  for (int j=0;j<J;j++){
    mu[j]=rmu[j];
  }
  sigma=*rsigma;
  normal=0; 
  priortype=*rpriortype;
  if (priortype==1) priortemp=rprior[0];
  if (*rmodel==2) normal=1;
   if (debug)  Rprintf("n%i J%i nseg%i level%i",n,J,nseg,level2); 
   if (level2==0) viterbi(data, lprob,  bestcp, rprior); 
   else viterbi_lev(data, lprob,  bestcp, bestlevel, rprior); 
}

void changepointci(double* cpconfint, double* cp, double* bestcp, int initsegci)
  { // obtain confidence interval of changepoint
    // start with cp with highest changepoint post prob, add adjacent points, whichever has higher post prob until total post prob within range is higher than entered ci

    double conf=0; // temporary cumulative confidence interval
    int ind=0,ind1,ind2; // ind: obs with highest prob of change-point, ind1: lower bound, ind2: higher bound
  
          if (verbose)    Rprintf("Estimating confidence intervals...\n");   
    for (int k=0;k<(nseg-1); k++) {
      if(initsegci){
          cpconfint[k]=seg[k+1];
          ind=seg[k+1]-1;
          if (debug)    Rprintf("index:%.5f\n",ind);   
      }
      else{  ind=bestcp[k];
             cpconfint[k]=ind;
          if (debug)    Rprintf("bestcp:%.5f\n",bestcp[k]);   
      }
      ind1=ind;
      ind2=ind;
      conf=cp[k*(n-1)+ind];
       
      while (conf<ci) {
        if ((ind1<=0)&(ind2>=(n-2))){ 
          if (debug)    Rprintf("a,"); 
           conf=1;  
           ind1=0;
           ind2=n-2;
           break;
	} else if (ind1<=0){
	  // increase upper bound 
          if (debug)    Rprintf("b,");  
	  ind2+=1;
	  conf+=cp[k*(n-1)+ind2];
	} else if (ind2>=n-2) {
	  // decrease lower bound 
          if (debug)    Rprintf("c,");  
	  ind1-=1;
	  conf+=cp[k*(n-1)+ind1];
	} else if ((cp[k*(n-1)+ind2+1]>cp[k*(n-1)+ind1-1])){
	  // increase upper bound 
          if (debug)    Rprintf("b,");  
	  ind2+=1;
	  conf+=cp[k*(n-1)+ind2];
	} else if ((cp[k*(n-1)+ind1-1]>cp[k*(n-1)+ind2+1])) {
	  // decrease lower bound 
          if (debug)    Rprintf("c,");  
	  ind1-=1;
	  conf+=cp[k*(n-1)+ind1];
	} else {
	  // change both if tied
          if (debug)    Rprintf("d,");  
          ind1-=1;
          ind2+=1;
          conf+=cp[k*(n-1)+ind1]+cp[k*(n-1)+ind2];
	}
      }
      cpconfint[nseg-1+k]=ind1+1;
      cpconfint[2*(nseg-1)+k]=ind2+1;
    }
  }


void postCPlev (double* data, int* rseg, int *nn, int *JJ, int *rnseg, double* lforward, double* lbackward, double* lprob, int *rprobs, double* cp, double* bestcp, double* bestlevel, double* rmu, double *rsigma, double* cpconfint, double* cpvector, double* levvector, double *rci, int *rviterbi, int *rinitsegci, int *rnsamples, int *rmodel, double* rprior, int *rpriortype, int *rmout, int *rlevelbased, int *rverbose, int *rdebug) {

int initsegci=*rinitsegci;
int doviterbi=*rviterbi;

verbose=*rverbose;
debug=*rdebug;
mout=*rmout;
n=*nn; 
J=*JJ;
nseg=*rnseg;

probs=*rprobs;
ci=*rci;
priortype=*rpriortype;
nsamples=*rnsamples;
normal=0;
if (priortype==1) priortemp=rprior[0];
if (*rmodel==2) normal=1;


if (debug)  Rprintf("normal? %i verbose %i\n",verbose);

 if (debug) Rprintf("n=%i,J=%i,seg=%i,\n",n,J,nseg);
seg.resize(nseg);
seg[0]=0;
for (int k=1;k<nseg;k++){
  seg[k]=rseg[k-1];
 if (debug) Rprintf("%i %i\n",k,seg[k]);
}
 mu.resize(J);

for (int j=0;j<J;j++){
  mu[j]=rmu[j];

}
sigma=*rsigma;
seg.push_back(n);

  // print the parameters
  if (debug) Rprintf("sigma_=%.8f\n",sigma);
  for (int j=0; j<J; j++)
    if (debug) Rprintf("mu[%i]=%.8f\n",j,mu[j]); 
  if (debug) Rprintf("\n"); 


  clock_t begin=clock();
  
  if (!mout) { 
    if (debug) Rprintf("creating matrices");
    double* lforward1=NULL;
    lforward1=new double[n*nseg*J];
    double* lbackward1=NULL;
    lbackward1=new double[n*nseg*J];
    if (debug) Rprintf("%i x %i x %i\n",n,nseg,J); 
    double* cp1=NULL;
    cp1=new double[(n-1)*(nseg-1)];
    if (debug) Rprintf("matrices created");
    forward_lev(data,lforward1,lprob,rprior);
    backward_lev(data,lbackward1,lforward1,lprob,cp1,rprior);

    if (doviterbi) viterbi_lev(data, lprob,  bestcp, bestlevel, rprior);
    if (ci>0) changepointci(cpconfint,cp1,bestcp,initsegci);
   if (nsamples>0) postCPsample(data, rmu, rsigma, lprob, rprobs, lbackward1, cpvector, levvector, rnsamples, nn, JJ, rnseg, rmodel,rprior,rpriortype,rlevelbased,rverbose,rdebug);
    delete [] lforward1;
    delete [] lbackward1;
    delete [] cp1;
    lforward1=NULL;
    lbackward1=NULL;
    cp1=NULL;
  }
  else {

    forward_lev(data,lforward,lprob,rprior);
    backward_lev(data,lbackward,lforward,lprob,cp,rprior);
    if (doviterbi)    viterbi_lev(data, lprob, bestcp, bestlevel,  rprior);
    if (ci>0) changepointci(cpconfint,cp,bestcp,initsegci);




    if (nsamples>0) postCPsample(data, rmu, rsigma, lprob, rprobs, lbackward, cpvector, levvector, rnsamples, nn, JJ, rnseg, rmodel,rprior,rpriortype,rlevelbased,rverbose,rdebug);
  }

   clock_t end=clock();
   if (debug)    Rprintf("Time elapsed: %.3f s\n", double(diffclock(end,begin)));

//  return 0;
}



void postCP (double* data, int* rseg, int *nn, int *JJ, double* lforward, double* lbackward, double* cp, double* bestcp, double* rmu, double *rsigma, double* lprob, int *rprobs, double* cpconfint, double* cpvector, double *rci, int *rviterbi, int *rinitsegci, int *rnsamples, int *rmodel, double* rprior, int *rpriortype, int *rmout, int *rlevelbased, int *rverbose, int *rdebug) {

int initsegci=*rinitsegci;
int doviterbi=*rviterbi;
verbose=*rverbose;
debug=*rdebug;
mout=*rmout;
n=*nn; 
J=*JJ;
nseg=J;
double* levvector=NULL;
levvector=new double[1];
levvector=0;
probs=*rprobs;
ci=*rci;
priortype=*rpriortype;
nsamples=*rnsamples;
normal=0;



if (priortype==1) priortemp=rprior[0];
if (*rmodel==2) normal=1;


if (debug)  Rprintf("normal? %i verbose %i\n",verbose);

 if (debug) Rprintf("n=%i,J=%i\n",n,J);
seg.resize(J);
seg[0]=0;
for (int j=1;j<J;j++){
  seg[j]=rseg[j-1];
 if (debug) Rprintf("%i %i\n",j,seg[j]);
}

seg.push_back(n);



  // print the parameters
 mu.resize(J);
for (int j=0;j<J;j++){
  mu[j]=rmu[j];
}
  sigma=*rsigma;
  if (debug) Rprintf("sigma=%.8f\n",sigma);
  for (int j=0; j<J; j++)
    if (debug) Rprintf("mu[%i]=%.8f\n",j,rmu[j]); 
  if (debug) Rprintf("\n"); 

  clock_t begin=clock();
  if (probs) initsegci=0;
  if (!mout) { 
    if (debug) Rprintf("creating matrices");
    double* lforward1=NULL;
    lforward1=new double[n*J];
    double* lbackward1=NULL;
    lbackward1=new double[n*J];
    if (debug) Rprintf("%i x %i\n",n,J); 
    double* cp1=NULL;
    cp1=new double[(n-1)*(J-1)];
    if (debug) Rprintf("matrices created");
    forward(data,lforward1,lprob,rprior);
    backward(data,lbackward1,lforward1,lprob,cp1,rprior);

    if (doviterbi) viterbi(data, lprob, bestcp, rprior);
    if (ci>0) changepointci(cpconfint,cp1,bestcp,initsegci);
    if (verbose) Rprintf("\n log(pevidence)=%.8f=%.8f\n",lforward1[0]+lbackward1[0],lforward1[n*J-1]+lbackward1[n*J-1]);
    if (nsamples>0) postCPsample(data, rmu, rsigma, lprob, rprobs, lbackward1, cpvector, levvector, rnsamples, nn, JJ, JJ, rmodel,rprior,rpriortype,rlevelbased,rverbose,rdebug);
    delete [] lforward1;
    delete [] lbackward1;
    delete [] cp1;
    lforward1=NULL;
    lbackward1=NULL;
    cp1=NULL;
  }
  else {
    forward(data,lforward,lprob,rprior);
    backward(data,lbackward,lforward,lprob,cp,rprior);
    if (doviterbi)    viterbi(data, lprob, bestcp, rprior);
    if (ci>0) changepointci(cpconfint,cp,bestcp,initsegci);

   if (verbose) Rprintf("\n log(pevidence)=%.8f=%.8f\n",lforward[0]+lbackward[0],lforward[n*J-1]+lbackward[n*J-1]);

    if (nsamples>0) postCPsample(data, rmu, rsigma, lprob, rprobs, lbackward, cpvector, levvector, rnsamples, nn, JJ, JJ, rmodel,rprior,rpriortype,rlevelbased,rverbose,rdebug);
  }

   clock_t end=clock();
   if (debug)    Rprintf("Time elapsed: %.3f s\n", double(diffclock(end,begin)));

//  return 0;
}





}


