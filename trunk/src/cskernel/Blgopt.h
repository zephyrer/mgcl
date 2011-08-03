/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// blgopt_ COMPUTES THE KNOTS T FOR THE OPTIMAL RECOVERY SCHEME OF ORDER  K 
//  FOR DATA AT  TAU(I), I=1,...,N . 
//
// ******  I N P U T  ****** 
//  TAU.....ARRAY OF LENGTH  N , CONTAINING THE INTERPOLATION POINTS. 
//    A S S U M E D  TO BE NONDECREASING, WITH TAU(I).LT.TAU(I+K),ALL I. 
//  N.....NUMBER OF DATA POINTS . 
//  K.....ORDER OF THE OPTIMAL RECOVERY SCHEME TO BE USED . 
//
// ******  W O R K  A R R A Y  ***** 
//  SCRTCH.....ARRAY OF LENGTH  (N-K)(2K+3) + 5K + 3 . THE VARIOUS 
//        CONTENTS ARE SPECIFIED IN THE TEXT BELOW . 
//
// ******  O U T P U T  ****** 
//  IFLAG.....INTEGER INDICATING SUCCESS (=1) OR FAILURE (=2) . 
//     IF IFLAG = 1, THEN 
//  T.....ARRAY OF LENGTH  N+K  CONTAINING THE OPTIMAL KNOTS READY FOR 
//        USE IN OPTIMAL RECOVERY. SPECIFICALLY,  T(1) = ... = T(K) = 
//        TAU(1)  AND  T(N+1) = ... = T(N+K) = TAU(N) , WHILE THE  N-K 
//        INTERIOR KNOTS  T(K+1), ..., T(N)  ARE CALCULATED AS DESCRIBED 
//        BELOW UNDER  *METHOD* . 
//     IF IFLAG = 2, THEN 
//        K .LT. 3, OR N .LT. K, OR A CERTAIN LINEAR SYSTEM WAS FOUND TO 
//        BE SINGULAR. 
//
// ******  M E T H O D  ****** 
//    THE (INTERIOR) KNOTS  T(K+1), ..., T(N)  ARE DETERMINED BY NEWTONS 
//  METHOD IN SUCH A WAY THAT THE SIGNUM FUNCTION WHICH CHANGES SIGN AT 
//   T(K+1), ..., T(N)  AND NOWHERE ELSE IN (TAU(1),TAU(N)) IS ORTHOGON- 
//  AL TO THE SPLINE SPACE  SPLINE( K , TAU )  ON THAT INTERVAL . 
//     LET  XI(J)  BE THE CURRENT GUESS FOR  T(K+J), J=1,...,N-K. THEN 
//  THE NEXT NEWTON ITERATE IS OF THE FORM 
//              XI(J)  +  (-)**(N-K-J)*X(J)  ,  J=1,...,N-K, 
//  WITH  X  THE SOLUTION OF THE LINEAR SYSTEM 
//                        C*X  =  D  . 
//  HERE,  C(I,J) = B(I)(XI(J)), ALL J, WITH  B(I)  THE I-TH B-SPLINE OF 
//  ORDER  K  FOR THE KNOT SEQUENCE  TAU , ALL I, AND  D  IS THE VECTOR 
//  GIVEN BY D(I) = SUM( -A(J) , J=I,...,N )*(TAU(I+K)-TAU(I))/K, ALL I, 
//  WITH  A(I) = SUM ( (-)**(N-K-J)*B(I,K+1,TAU)(XI(J)) , J=1,...,N-K ) 
//  FOR I=1,...,N-1, AND  A(N) = -.5 . 
//    (SEE CHAPTER  XIII  OF TEXT AND REFERENCES THERE FOR A DERIVATION) 
//    THE FIRST GUESS FOR  T(K+J)  IS  (TAU(J+1)+...+TAU(J+K-1))/(K-1) . 
//     ITERATION TERMINATES IF  MAX(ABS(X(J))) .LT. T O L  , WITH 
//                 T O L  =  T O L R T E *(TAU(N)-TAU(1))/(N-K) , 
//  OR ELSE AFTER  N E W T M X  ITERATIONS , CURRENTLY, 
//                 NEWTMX, TOLRTE / 10, .000001 
//
//     DIMENSION SCRTCH((N-K)*(2*K+3)+5*K+3), T(N+K) 
// URRENT FORTRAN STANDARD MAKES IT IMPOSSIBLE TO SPECIFY THE PRECISE 
// DIMENSIONS OF SCRTCH AND T WITHOUT THE INTRODUCTION OF OTHERWISE 
//  SUPERFLUOUS ADDITIONAL ARGUMENTS . 
void blgopt_(const double *tau, int n, int k, double *scrtch, double *t, int *iflag);
