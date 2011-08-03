/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
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
void b2vb_(int k, int npk,const double *t, double x, int left, double *biatx);
