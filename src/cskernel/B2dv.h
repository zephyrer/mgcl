/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
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
//                 ZERO DIVISION OCCURS. 
//    JDERIV       INDICATES ORDER OF DERIVATIVE , MUST BE NON-NEGATIVE 
//                 , MAY BE ZERO. 
// ***OUTPUT* 
//    RDATX(K)     EVALUATED COEFFICIENTS OF SUPPOSEDLY NON ZERO 
// ***NOTE* 
//    RDATX(I), 1<=I<=K ARE THE COEFFICIENTS OF RCOEF(J) 
//             LEFT-K+1 <= J <= LEFT , 
//    I.E. THE COEFFICIENT OF RCOEF(LEFT-K+I) IS RDATX(I). 
void b2dv_(int k, int npk, const double *t, double x, int left, int jderiv, double *rdatx);
