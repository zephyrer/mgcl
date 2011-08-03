/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <malloc.h>
#include "cskernel/b2dv.h"
#include "cskernel/bk2fli.h"

// REAL FUNCTION TO EVALUATE LEFT-CONTINUOUS JDERIV-TH DERIVATIVE 
// AT THE PARAMETER VALUE X OF THE B-REP. (T,RCOEF). 
// *** INPUT * 
//     K,N,T(N+K),RCOEF(N).......B-REP TO EVALUATE,  ORDER, B-REP DIM- 
//              EMNSION, KNOT VECTOR, AND B-COEFFICIENTS EACH. 
//     X        VALUE AT WHICH THE DERIVATIVE IS EVALUATED. 
//     JDERIV   ORDER OF THE DERIVATIVE,MAY BE ZERO. 
// *** OUTPUT * 
//     BLEL   THE VALUE OF THE B-SPLINE AT THE PARMETER X. 
// *** NOTE * 
//     FUNCTION BLER EVALUATES RIGHT-CONTINUOUS DERIVATIVE WHILE 
//     BLEL DOES LEFT-CONTINUOUS. 
double blel_(int k, int n, const double *t, const double *rcoef, double x, int jderiv){
    double ret_val;
    int i;
    int id, ki, npk;

    double rdatx_buf[20];
	double* rdatx=rdatx_buf;
	if(k>20)
		rdatx=(double*)malloc(sizeof(double)*(k));

    // Parameter adjustments 
    --rcoef;
    --t;

    // Function Body 
    npk = n+k;
    ki=bk2fli_(npk, &t[1], x);
    if(x>t[1] && x<=t[npk]){
		while(x==t[ki])
			--ki;
    }
    b2dv_(k, npk, &t[1], x, ki, jderiv, rdatx);
    ret_val = 0.f;
    ki -= k;
    for(i=1; i<=k; ++i){
		id = ki+i;
		if(id<1) id=1;
		if(id>n) id=n;
		ret_val += rcoef[id]*rdatx[i-1];
    }

	if(k>20) free(rdatx);
    return ret_val;
}
