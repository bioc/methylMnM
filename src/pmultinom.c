# include <R.h>
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <Rdefines.h>
#include <Rinternals.h>
void pcalculate(double *T, int *SIZE, int *length, double *P1, double *P2, double *P3, double *P4, double *C1, double *C2, double *pvalue)
{
   int ii;

   for(ii=0; ii<*length; ii++)
   {
       //Multinom(*T[ii], *SIZE[ii], *P1[ii], *P2[ii], *P3[ii], *P4[ii], *C1, *C2);          
       //pMultinom(double t, int size, double p1, double p2, double p3, double p4, double c1, double c2)  
	   double t=T[ii], p1=P1[ii], p2=P2[ii], p3=P3[ii], p4=P4[ii], c1=*C1, c2=*C2;
       int size=SIZE[ii];

        int i, j, k, r; 
       	int i1,j1,k1,r1;
     	double p=0, pp, t1;
      for(i=0; i<=size; i++){
       for(j=0; j<=size-i; j++){ 
          for(r=0; r<=size-i-j; r++){
             double sj1=1, sr1=1, sk1=1, si1=1;
              k=size-i-j-r;
			  t1=c1*i*j-c2*r*k;
			  if (t1<0) t1=-t1;
              if(t1>=t) 
			  { if(i==0){si1=1;} 
			    else{
				 for(i1=1; i1<=i; i1++)
				 { 
					si1=p1*si1;
				 }
				}
	
			if(j==0){sj1=1;} 
			else{
			   for(j1=i+1; j1<=i+j; j1++)
			   { sj1=sj1*(j1*p4/(j1-i));
			   }
			}
          	if(r==0){sr1=1;} 
			else{
               for(r1=i+j+1; r1<=i+j+r; r1++)
				{ sr1=sr1*(r1*p2/(r1-i-j));
			   }
			}
			if(k==0){sk1=1;} 
			else{
			   for(k1=i+j+r+1; k1<=size; k1++)
				 { sk1=sk1*(k1*p3/(k1-i-j-r));
			   }
			}
                
				pp=si1*sj1*sr1*sk1;
				p=p+pp;
	
			       }
			else
			{p=p+0;}
		  }
	   }
   }
         pvalue[ii]=p;	   
	   if(fmod(ii,50000)==0) printf("%d\n",ii);
   }
}



