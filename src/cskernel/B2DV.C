/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/b2vb.h"

//     B2DV EVALUATES THE K COEFFICIENTS OF B-COEFFICIENTS ,GIVEN 
//     THEIR KNOT VECTORS. 
// ***INPUT* 
//    K,NPK,T(NPK)......PROVIDE KNOT VECTOR OF ORDER K AND LENGTH NPK. 
//    X            PARAMETER VALUE AT WHICH DERIVATIVE TO EVALUATE 
//    LEFT         SPECIFIES WHERE X IS LOCATED IN T(.), I.E. 
//                 T(LEFT)<= X <T(LEFT+1) IN GENERAL. 
//                 HOWEVER, WHEN X IS OUTSIDE THE RANGE OF T(.), THIS 
//                 IS NOT THE CASE. 
//     CASE 1) WHEN X IS GREATER THAN T(NPK), T(LEFT) < T(LEFT+1) < X. 
//     CASE 2) WHEN X IS LESS THAN T(1), X < T(LEFT) < T(LEFT+1). 
//                 T(LEFT+1) MUST BE GREATER THAN T(LEFT), IF EQUAL 
//                 ZERO DIVISION OCCURS. 
//    JDERIV       INDICATES ORDER OF DERIVATIVE , MUST BE NON-NEGATIVE 
//                 , MAY BE ZERO. 
// ***OUTPUT* 
//    RDATX(K)     EVALUATED COEFFICIENTS OF SUPPOSEDLY NON ZERO 
// ***NOTE* 
//    RDATX(I), 1<=I<=K ARE THE COEFFICIENTS OF RCOEF(J) 
//             LEFT-K+1 <= J <= LEFT , 
//    I.E. THE COEFFICIENT OF RCOEF(LEFT-K+I) IS RDATX(I). 
void b2dv_(int k, int npk, const double *t, double x, int left, int jderiv, double *rdatx){
    int kmjd, i, j;
    double term, saved,fj;
    int i1, i2, jp1, lmj;
// *******************START OF B2DV*********************************** 
    // Parameter adjustments 
    --rdatx;
    --t;

    // Function Body 
    for (i=1; i<=k; ++i)
		rdatx[i]=0.f;
    if (jderiv>=k)
		return;

    kmjd = k-jderiv;
    b2vb_(kmjd, npk, &t[1], x, left, &rdatx[1]);
// NOW B-SPLINE RDATX(I) 1<=I<=K-JDERIV(ORDER OF K-JDERIV) ARE OBTAINED 
// GET DERIVATIVE PART OF THE COEFFICIENTS 
    j = kmjd;
	while(j<k){
	    jp1 = j+1;
		saved = 0.f;
	    fj = (double) j;
		lmj = left-j;
		for(i=1; i<=j; ++i){
			i1=lmj+i;
			if(i1<1)
				i1=1;
			i2=left+i;
			if(i2>npk)
				i2=npk;
			term=rdatx[i]*fj / (t[i2]-t[i1]);
			rdatx[i] = saved-term;
			saved =term;
	    }
		rdatx[jp1] = saved;
	    j = jp1;
	}
}
