# include <R.h>
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <Rdefines.h>
#include <Rinternals.h>
void cpgcount(int *data2, int *data3, int *cpg2, int *cpg3, int *datalength, int *cpglength,  int *count)
{
   int i, j, k=0,count2=0;

      for(i=0; i<*cpglength; i++){
		  int count1=0;
		  k=k-count2;
       for(j=0; j<*datalength; j++){ 
		   if (data2[j]<=cpg2[i] && data3[j]>=cpg3[i]) 
			   {count1=count1+1;		   
		   }
		   else{
			   if(data2[j]>=cpg2[i])
			   {
				   k=j;
				   break;
			   } 
		   }
	   }
	   count2=count1;
          count[i]=count1;
	   if(fmod(i,200000)==0) printf("%d\n",i);	  
	  }
}



