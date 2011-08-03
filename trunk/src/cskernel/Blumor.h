/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLUMOR         BLUMOR FOR BLUMOV 
// BLUMOR GENERATES RATIO'S OF B-COEFFICIENT TRANSLATION. 
// *** INPUT * 
//       KMOVI.....INDICATES WHAT KIND OF MOVE IS BEING PERFORMED. 
//               =1 : TWO POINTS FIXED AND CENTER OF MOVE BETWEEN THEM. 
//               =2 : ONE POINT,T(I),FIXED. THE OTHER POINT IS FREE 
//               =3 : ONE POINT,T(J+1),FIXED. THE OTHER POINT IS FREE 
//       I,J....GIVE ID OF KNOT VECTOR T(.) THAT SHOULD BE FIXED. 
//                 T(I)<= TAU <=T(J+1)   K<=I<=N, K<=J<=N. 
//       TAU......IS THE CENTER OF TRANSLATION, PARAMETER VALUE OF THE 
//                B-REP. 
//       K,N,T....DESCRIBE KNOT VECTOR OF THE B-REP. 
//                K : ORDER,      N : B-REP DIMENSION 
//                T(N+K) : KNOT VECTOR 
// *** OUTPUT * 
//       I,J.....UPDATED I AND J. 
//       ISTRT,IEND......START AND END ID OF B-COEF RATIO(.,.) THAT 
//                SHOULD BE MODIFIED. 
//       RATIO(M)...RATIO(M) CONTAINS RATIO OF TRANSLATION 
//                  RATIO(M) IS FOR RCOEF(M,.) 
//                    ISTRT<= M <= IEND 
// ***WORK* 
//       RATIO(N)...WORK AREA OF LENGTH N. 
void blumor_(int kmovi, int *i, int *j, 
	double tau, int k, int n, const double *t, int *istrt,
	int *iend, double *ratio);

