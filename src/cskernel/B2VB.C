/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <malloc.h>

#if defined(_DEBUG)
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//The original codes of this program comes from the FORTRAN code BSPLVB of
//"A Practical Guide to Splines" by Carl de Boor.

// CALCULATES THE VALUE OF ALL POSSIBLY NONZERO B-SPLINES AT  X 
//  OF ORDEER  K  WITH KNOT SEQUENCE T(.). 
// ******  I N P U T  ****** 
//   K,NPK,T(NPK)....PROVIDE KNOT VECTOR OF ORDER K AND 
//                    KNOT LENGTH NPK. 
//  X.....THE POINT AT WHICH THE B-SPLINES ARE TO BE EVALUATED. 
//  LEFT.....AN INTEGER WHICH HOLDS THE FOLLOWING CONDITIONS; 
//           T(1) <= X <= T(NPK), T(LEFT) <= X <= T(LEFT+1) 
//                         AND    T(LEFT) < T(LEFT+1). 
//             1 <= LEFT <= NPK-1. 
// .....W A R N I N G ..... WHEN T(LEFT)=T(LEFT+1), ZERO DIVISION CAUSED 
// ..... THESE CONDITIONS HOLD AS LONG AS BK2FLI IS USED TO FIND LEFT, 
//        SEE BK2FLI FOR DETAIL. 
// ******  O U T P U T  ****** 
//  BIATX(K).....ARRAY OF LENGTH K, WITH  BIATX(I) 
//        CONTAINING THE VALUE AT  X  OF THE POLYNNMIAL OF ORDER 
//        K  WHICH AGREES WITH THE B-SPHINE  B(LEFT-K+I,K,T) 
//        ON THE INTERVAL (T(LEFT),T(LEFT+1)) . 
// ******  M E T H O D  ****** 
//  THE RECURRENCE RELATION 
//                       X - T(I)              T(I+J+1) - X 
//     B(I,J+1)(X)  =  -----------B(I,J)(X) + ---------------B(I+1,J)(X) 
//                     T(I+J)-T(I)            T(I+J+1)-T(I+1) 
//  IS USED (REPEATEDLY) TO GENERATE THE (J+1)-VECTOR 
//  B(LEFT-J,J+1)(X),...,B(LEFT,J+1)(X)  FROM THE J-VECTOR 
//  B(LEFT-J+1,J)(X),&..,B(LEFT,J)(X), STORING THE NEW VALUES IN 
//  BIATX  OVER THE OLD. THE FACTS THAT 
//            B(I,1) = 1  IF  T(I) .LE. X .LT. T(I+1) 
//  AND THAT 
//            B(I,J)(X) = 0  UNLESS  T(I) .LE. X .LT. T(I+J) 
//  ARE USED. 
// ****** N O T E ****** 
// 1) WHEN T(1) <= X  AND LEFT < K, T(1) IS USED AS MISSING KNOT(S) AS 
//    MANY AS NECESSARY, I.E.  T(1),.....,T(1),T(2),....   AS KNOTS. 
// 2) WHEN T(NPK) >=  X  AND LEFT > NPK-K, T(N+K) IS USED AS MISSING 
//    KNOT(S) AS MANY AS NECESSARY, I.E. 
//             T(NPK-1),T(NPK),....,T(NPK) AS NEW KNOT CONFIG. 
void b2vb_(int k, int npk,const double *t, double x, int left, double *biatx){
    double term;
    int i, j;
    double saved;
    int ll, lr;
    int jp1;

    double deltal_buf[20], deltar_buf[20];
	double* deltal=deltal_buf; double* deltar=deltar_buf;
	if(k>20){
		deltal=(double*)malloc(sizeof(double)*k);
		deltar=(double*)malloc(sizeof(double)*k);
	};
    // Parameter adjustments 
    --biatx;
    --t;

    // Function Body 
    j = 1;
    biatx[1] = 1.f;
	while(j < k) {
	    jp1 = j + 1;
		lr = left + j;
		if (lr > npk) lr = npk;
		deltar[j-1] = t[lr] - x;
	    ll = left+1-j;
		if (ll < 1) ll = 1;
		deltal[j-1] = x-t[ll];
		saved = 0.f;
		for	(i = 1; i <= j; ++i) {
			term = biatx[i]/(deltar[i-1] + deltal[jp1-i-1]);
			biatx[i] = saved+deltar[i-1]*term;
			saved = deltal[jp1-i-1] * term;
		}
		biatx[jp1] = saved;
		j = jp1;
	}

	if(k>20){
		free(deltal); free(deltar);
	}
}
