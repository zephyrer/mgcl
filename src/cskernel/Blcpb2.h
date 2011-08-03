/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLCPB2 COMPUTES B-COEFFICIENTS OF B-REP, FROM P-REP.(RBRK,PCOEF,L) AND KNOT VECTOR T. 
// ***INPUT* 
//    PCOEF(K,IPC,NCD).PP-COEFFICIENTS , I.E. 
//             PCOEF(J,I,,)=D**(J-1)(F(RBRK(I))   1<=I<=L. 
//    L................INDICATES NUM OF INTERVAL OF PP-REP. 
//    K................ORDER OF PP-REP(B-REP) 
//    N................B-REP DIMENSION 
//    T(N+K)...........KNOT VECTOR OF B-REP 
//    NCD..............SPACE DIMENSION OF THE PP-REP 
//    IPC..............LENGTH OF 2ND DIMENSIONED ARRAY OF THE VARIALBLE 
//                     PCOEF, MUST BE .GE.L+1 . 
//    IRC..............ROW DIMENSION OF THE VARIABLE RCOEF 
// ***OUTPUT*** 
//    RCOEF(IRC,NCD)...B-COEFFICIENTS 
void blcpb2_(const double *pcoef, int l, int k, int n,
			const double *t, int ncd, int ipc, int irc, double *rcoef);
