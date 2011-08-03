/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLUDKT GENERATES B-COEFFICIENTS OF DIFFERENT KNOT CONFIGURATION T2. 
// *** INPUT * 
//     K,N1,T1(N1+K),RCOE1(IRC1,NCD)IRC1,NCD...DESCRIBE THE OLD B-REP 
//          T1:KNOT VECTOR , RCOE1:B-COEF , N1:B-REP DIMENSION 
//     N2,T2(N2+K)....NEW BREP DIMENSION AND KNOT VECTOR 
//     IRC2....ROW DIMENSION OF THE VARIABLE RCOE2 
// *** OUTPUT * 
//     RCOE2(IRC2,NCD)..THE NEW B-COEF OBTAINED 
//     IFLAG... =1 :SUCCESFUL RETURN, <>1 : FAILURE BECAUSE OF ILLEGAL 
//                                          KNOT VECTOR T2. 
// *** WORK * 
//     WK1(N2),WK2(N2,9) 
// *** NOTE * 
//     THE NEW KNOT T2 MUST SATISFY THE FOLLOWLING CONDITION; 
//          T1(K)<=T2(K), AND T2(N2+1)<=T1(N1+1) 
void bludkt_(int k, int n1, const double *t1, 
	const double *rcoe1, int irc1, int ncd, int n2, 
	const double *t2, int irc2, double *wk1, double *wk2, 
	double *rcoe2, int *iflag);
