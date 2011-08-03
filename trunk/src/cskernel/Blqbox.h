/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//BLQBOX computes box of given line B-rep, i.e. box of B-coefficients of the
// B-Rep. 
// Given TS(start parameter) and TE(end parameter), computes box boundary 
// RBOUND of the partial B-rep. 
// *** INPUT  * 
//   K,N,T(N+K),RCOEF(IRC,NCD),IRC.......PROVIDE B-REP OF 
//               ORDER K,B-REP DIMENSION N, KNOT VECTOR T(.), AND 
//               B-COEFF'S RCOEF(.,.). 
//   NCD........SPACE DIMENSION OF THE B-REP. 
//   TS,TE......indicates the paramter range o the partal B-REP. 
// *** OUTPUT * 
//   RBOUND(2,NCD).....the box boundary, RBOUND(1,.): Minimum Point and 
//                                       RBOUND(2,.): Maximum Point. 
// *** WORK   * 
//     TW(N+K),BATJ(K,K)..... OF EACH LENGTH 
void blqbox_(int kp, int n, const double *t, 
	const double *rcoef, int irc, int ncd, double ts, 
	double te, double *tw, double *batj, double *rbound);
