#include <R.h>
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <Rdefines.h>
#include <Rinternals.h>
void pvalueclassify(int *type1, int *type2, int *type3, int *type4, int *sm1chring1, int *sm1chring2, int *sm1chring3, int *sm1chring4,double *p, int *typelength, int *sm1chringlength,  double *pvalue)
{
   int i,j;  
      for(i=0; i<*sm1chringlength; i++){
        for(j=0; j<*typelength; j++ ){ 
		   if (sm1chring1[i]==type1[j] && sm1chring2[i]==type2[j] &&sm1chring3[i]==type3[j] && sm1chring4[i]==type4[j]) {  
			   pvalue[i]=p[j];	
		       break;
		   }
		}
	    if(fmod(i,200000)==0) printf("%d\n",i);
	  }
}
