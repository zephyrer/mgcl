/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// CONVERTS THE B-REPRESENTATION  K,N,T,BCOEF OF SOME SPLINE INTO ITS 
//  PP-REPRESENTATION  tau, COEF, L, K . 
// ******  I N P U T  ****** 
//  ORDER.....ORDER OF THE SPLINE 
//  N.....LENGTH OF  BCOEF  AND  DIMENSION OF SPLINE SPACE  SPLINE(K,T) 
//  T.....KNOT SEQUENCE, OF LENGTH  N+K 
//  BCOEF(irc,ncd).....B-SPLINE COEFFICIENT SEQUENCE, OF LENGTH  N,
//				and space dimension ncd
//The row dimension of BCOEF is irc.
//ipc.....indicates the structure of coef as described in coef.
// ******  W O R K   A R E A  ****** 
//  SCRTCH......OF SIZE  (K,K,ncd) , NEEDED TO CONTAIN BCOEFFS OF A PIECE OF 
//                THE SPLINE AND ITS  K-1  DERIVATIVES 
// ******  O U T P U T  ****** 
//  TAU.....BREAKPOINT SEQUENCE, OF LENGTH  L+1, CONTAINS (IN INCREAS- 
//        ING ORDER) THE DISTINCT POINTS IN THE SEQUENCE  T(K),...,T(N+1) 
//  COEF(K,ipc,ncd).....ARRAY OF SIZE (K,ipc,ncd), WITH  COEF(I,J, m) =
//		(I-1)ST DERIVATIVE of m-th space dimension OF SPLINE AT tau(J)
//		FROM THE RIGHT                                                  
//  L...NUMBER OF POLYNOMIAL PIECES WHICH MAKE UP THE SPLINE IN THE
//		INTERVAL  (T(K), T(N+1)). L's maximum is N+1-K.                 
// ******  M E T H O D  ****** 
//   FOR EACH BREAKPOINT INTERVAL, THE  K  RELEVANT B-COEFFS OF THE 
//   SPLINE ARE FOUND AND THEN DIFFERENCED REPEATEDLY TO GET THE B-COEFFS 
//   OF ALL THE DERIVATIVES OF THE SPLINE ON THAT INTERVAL. THE SPLINE AND 
//   ITS FIRST  K-1  DERIVATIVES ARE THEN EVALUATED AT THE LEFT END POINT 
//   OF THAT INTERVAL. 
void blcbpn_(
	int order, int n, const double *t,const double *bcoef, int irc, int ncd, int ipc,
	double *scrtch,	double *tau, double *coef, int *l);
