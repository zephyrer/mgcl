/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// CONVERTS THE B-REPRESENTATION  K,N,T,BCOEF OF SOME SPLINE INTO ITS 
//  PP-REPRESENTATION  BREAK, COEF, L, K . 
// ******  I N P U T  ****** 
//  K.....ORDER OF THE SPLINE 
//  N.....LENGTH OF  BCOEF  AND  DIMENSION OF SPLINE SPACE  SPLINE(K,T) 
//  T.....KNOT SEQUENCE, OF LENGTH  N+K 
//  BCOEF.....B-SPLINE COEFFICIENT SEQUENCE, OF LENGTH  N 
// ******  W O R K   A R E A  ****** 
//  SCRTCH......OF SIZE  (K,K) , NEEDED TO CONTAIN BCOEFFS OF A PIECE OF 
//                THE SPLINE AND ITS  K-1  DERIVATIVES 
// ******  O U T P U T  ****** 
//  BREAK.....BREAKPOINT SEQUENCE, OF LENGTH  L+1, CONTAINS (IN INCREAS- 
//        ING ORDER) THE DISTINCT POINTS IN THE SEQUENCE  T(K),...,T(N+1) 
//  COEF.....ARRAY OF SIZE (K,N), WITH  COEF(I,J) = (I-1)ST DERIVATIVE 
//   OF  SPLINE AT BREAK(J) FROM THE RIGHT 
//  L...NUMBER OF POLYNOMIAL PIECES WHICH MAKE UP THE SPLINE IN THE IN- 
//        TERVAL  (T(K), T(N+1)) 
// ******  M E T H O D  ****** 
//     FOR EACH BREAKPOINT INTERVAL, THE  K  RELEVANT B-COEFFS OF THE 
//   SPLINE ARE FOUND AND THEN DIFFERENCED REPEATEDLY TO GET THE B-COEFFS 
//   OF ALL THE DERIVATIVES OF THE SPLINE ON THAT INTERVAL. THE SPLINE AND 
//   ITS FIRST  K-1  DERIVATIVES ARE THEN EVALUATED AT THE LEFT END POINT 
//   OF THAT INTERVAL. 
void blcbp_(
	int k, int n, const double *t, const double *bcoef,
	double *scrtch, double *tau, double *coef, int *l
);
