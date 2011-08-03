/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <malloc.h>
#include "cskernel/b2dv.h"
#include "cskernel/bk2fli.h"

//   B1A COMPUTES THE INTEGRAL OF B-REP (K,N,T,RCOEF) FROM PARAMETER 
//   T(1) TO X. 
//***INPUT*
//  K,N,T(N+K),RCOEF(N)......PROVIDE B-REP(OF ONE SPACE DIMENSION), 
//             FUNCTION TO INTEGRATE. 
//  X.....SPECIFIY PARAMETER RANGE OF THE INTEGRATION, FROM T(1) TO X.
//      X MUST LIE IN THE PARAMETER RANGE OF THE B-REP, I.E. 
//                    T(1) <= X <= T(N+K). 
//***OUTPUT*
//  B1A......INTEGRAL FROM PARAMETER T(1) TO X. 
double b1a_(int k, int n,const double *t,const double *rcoef, double x){
    // System generated locals 
    double ret_val;

    // Local variables 
    int jxmk, i, j;
    int jc;
    double fk;
    int jt1, jt2, jx;
    double bci;
    int npk;

// *******************START OF B1A**********************************
    double rdatx_buf[21];
	double* rdatx=rdatx_buf;
	// type cast (void* -> double*)
	if(k>20) rdatx=(double*)malloc(sizeof(double)*k+1);

    // Parameter adjustments 
    --rcoef;
    --t;

    // Function Body 
    npk = n+k;
    if(x<t[1]) x=t[1];
    if(x>t[npk]) x=t[npk];
    jx=bk2fli_(npk, &t[1], x);
    b2dv_(k+1, npk, &t[1], x, jx, 0, rdatx);
    jxmk = jx-k;
    fk = 1.f/(double)k;
    bci = 0.f;
    for(j=2-k; j<=jxmk; ++j){
		jt1 = j;
		if(jt1<1) jt1=1;
		bci+=rcoef[jt1]*(t[j+k]-t[jt1]);
    }
    ret_val = bci*rdatx[0]*fk;
    for(i=1; i<=k; ++i){
		j = jxmk + i;
		jt1 = j;
		if(jt1<1) jt1=1;
		jt2 = j+k;
		if(jt2>npk) jt2=npk;
		jc=jt1;
		if(jc>n) jc=n;
		bci += rcoef[jc]*(t[jt2]-t[jt1]);
		ret_val += bci*rdatx[i]*fk;
	}

	if(k>20) free(rdatx);
    return ret_val;
}
