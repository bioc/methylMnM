#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "methylmnm.h"
         
//Register C routines
R_CMethodDef cMethods[]={
	{"calculatecount", (DL_FUNC)&calculatecount, 7},
	{"calculatecount1", (DL_FUNC)&calculatecount1, 7},
	{"calculatecountneg", (DL_FUNC)&calculatecountneg, 7},
	{"cpgcount", (DL_FUNC)&cpgcount, 7},
	{"pcalculate", (DL_FUNC)&pcalculate, 10},
    {"pvalueclassify", (DL_FUNC)&pvalueclassify, 12},
	{NULL, NULL, 0}
};

void R_init_methylMnM(DllInfo *info)
{
	R_registerRoutines(info, cMethods, NULL, NULL, NULL);
}
