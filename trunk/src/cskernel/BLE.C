/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <malloc.h>
#include "cskernel/b2dv.h"
#include "cskernel/bk2fli.h"

// BLE EVALUATES THE JDERIV-TH DERIVATIVE(S), GIVEN B-COEF'S 
// RCOEF IN COLUMN-WISE OR ROW-WISE. B-COEFFICIENTS RCOEF MAY 
// BE MULTIPLE OF THE SAME KNOT CONFIGURATION. 
// ***INPUT* 
//    K,N,T(N+K),RCOEF(IRC,.),IRC,NCD....DESCRIBE B-REP. 
//           ORDER, B-REP DIMENSION, KNOT VECTOR, B-COEF, AND SPACE 
//           DIMENSION. 
//             RCOEF(I,J)   1<=I<=N  (J=1...NCD) 
//    X,JDERIV....PARAMETER VALUE AT WHICH DERIVATIVE TO EVALUATE AND 
//           ORDER OF DERIVATIVE. THE DERIVATIVE ORDER JDERIV MUST BE 
//           NON-NEGATIVE, MAY BE ZERO. 
// ***OUTPUT* 
//    P(NCD)         EVALUATED DERIVATIVE(S) 
// ***NOTE* 
//    BLE IS EASY-TO-USE VERSION OF BLEVAL, AND HAVE BETTER PERFORMANCE 
//    THAN BLEVAL. 
void ble_(
	int k,int n,const double *t,const double *rcoef,int irc,int ncd,double x,int jderiv,double *p
){
    // Local variables 
    int irmk, i, l;
    int id;
    int irm, npk;

// *******************START OF BLEVAL***********************************
    double rbatir_buf[20];
	double* rbatir=rbatir_buf;
	if(k>20)
		rbatir=(double*)malloc(sizeof(double)*k);

    // Parameter adjustments 
    --p;
    rcoef -= irc+1;

    // Function Body 
    npk = n+k;
//    FIND WHERE X IS LOCATED IN T(.). 
    irm=bk2fli_(npk, t, x);
    b2dv_(k, npk, t, x, irm, jderiv, rbatir);
    irmk = irm-k;
//    NOW COEFFICIENTS OF B-COEF'S ARE OBTAINED, 
//    GET JDERIV-TH DERIVATIVES BY MULTIPLYING THEM. 
    for(l = 1; l <= ncd; ++l)
		p[l] = 0.f;
    for(i = 1; i <= k; ++i){
		id = irmk+i;
		if(id<1)
			id=1;
		if(id>n)
			id=n;
		for(l=1; l<=ncd; ++l)
			p[l] += rcoef[id+l*irc]*rbatir[i-1];
    }

	if(k>20)
		free(rbatir);
}
