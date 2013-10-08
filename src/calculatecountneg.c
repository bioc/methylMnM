# include <R.h>
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <Rdefines.h>
#include <Rinternals.h>
void calculatecountneg(int *data2, int *data3, int *cpg2, int *cpg3, int *datalength, int *cpglength,  int *count)
{
   int i,j, k=0,count2=0;  
      for(i=0; i<*cpglength; i++){
		  int count1=0;
          k=k-count2;
       for(j=k; j<*datalength; j++ ){ 
		   if (data3[j]<=cpg3[i] && data3[j]>cpg2[i]) 
		   {count1=count1+1;		   
		   }
		   else{
			   if(data3[j]>cpg3[i])
			   {
				   k=j;
				   break;
			   } 
		   }

	   }
      // k=k+count1+count2;
	   count2=count1;
       count[i]=count1;
	   if(fmod(i,100000)==0) printf("%d\n",i);
   }
}



