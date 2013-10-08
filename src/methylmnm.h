#include <string.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

/* Functions called by R code */

void calculatecount1(int *data2, int *data3, int *cpg2, int *cpg3, int *datalength, int *cpglength,  int *count);
void calculatecount(int *data2, int *data3, int *cpg2, int *cpg3, int *datalength, int *cpglength,  int *count);
void calculatecountneg(int *data2, int *data3, int *cpg2, int *cpg3, int *datalength, int *cpglength,  int *count);
void cpgcount(int *data2, int *data3, int *cpg2, int *cpg3, int *datalength, int *cpglength,  int *count);
void pcalculate(double *T, int *SIZE, int *length, double *P1, double *P2, double *P3, double *P4, double *C1, double *C2, double *pvalue);
void pvalueclassify(int *type1, int *type2, int *type3, int *type4, int *sm1chring1, int *sm1chring2, int *sm1chring3, int *sm1chring4,double *p, int *typelength, int *sm1chringlength,  double *pvalue);
