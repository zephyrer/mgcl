/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <malloc.h>
#include "cskernel/b2dv.h"
#include "cskernel/bk2fli.h"

//   REAL FUNCTION TO EVALUATE JDERIV-TH DERIVATIVE AT THE PARAMETER 
//   VALUE X OF THE B-REP (N,T,RCOEF). 
// *** INPUT *
//     K,N,T(N+K),RCOEF(N)......B-REP TO EVALUATE,    ORDER, B-REP DIM- 
//     X        VALUE AT WHICH THE DERIVATIVE IS EVALUATED. 
//     JDERIV   ORDER OF THE DERIVATIVE,MAY BE ZERO. 
// *** OUTPUT *
//     BLER    THE VALUE OF THE B-SPLINE AT THE PARMETER X. 
double bler_(int k, int n, const double *t, const double *rcoef, double x, int jderiv){
    int i, id, ki, npk, kimk;
    double ret_val;

    double rdatx_buf[20];
	double* rdatx=rdatx_buf;
	if(k>20)
		rdatx=(double*)malloc(sizeof(double)*(k));

    // Parameter adjustments 
    --rcoef;

    npk = n+k;
	// ***     FKIND KI S.T.  T(KI) <= X < T(KI+1). 
    ki=bk2fli_(npk, t, x);
    b2dv_(k, npk, t, x, ki, jderiv, rdatx);
    kimk = ki-k;
    ret_val = 0.f;
    for(i=1; i<=k; ++i){
		id = kimk + i;
		if(id<1)
			id = 1;
		if(id>n){
		    id=n;
		}
		ret_val += rcoef[id]*rdatx[i-1];
    }

	if(k>20) free(rdatx);
    return ret_val;
}
