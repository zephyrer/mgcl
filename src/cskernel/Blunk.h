/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLUNK GENERATES B-COEFFICIENTS OF NEW KNOT CONFIGURATION T2. 
// *** INPUT * 
//     K,N1,T1(N1+K),RCOE1(IRC1,NCD),IRC1,NCD...DESCRIBE THE OLD B-REP 
//     N2,T2(N2+K).....KNOT VECTOR OF THE NEW B-REP 
//     IRC2....ROW DIMENSION OF THE VARIABLE RCOE2 
// *** OUTPUT * 
//     RCOE2(IRC2,NCD)..THE NEW B-COEF OBTAINED 
// *** WORK * 
//     BATJ(K,K)  LENGTH OF K*K 
// *** NOTE * 
//     THE NEW KNOT T2 MUST SATISFY THE FOLLOWLING CONDITION; 
//          FOR ANY I (1<=I<=N2-1), THE NUM OF J S.T. 
//          T2(I)<= T1(J) < T2(I+1) 
//          MUST BE LESS OR EQUAL 1. 
void blunk_(int k, int n1, const double *t1, 
	const double *rcoe1, int irc1, int ncd, int n2, 
	const double *t2, int irc2, double *batj, double *rcoe2);

