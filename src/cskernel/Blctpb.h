/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLCTPB COMOUTES B-COEFFICIENTS OF THE B-REP OF PP-REP(RBRK,PCOEF,L). 
// ***INPUT* 
//    RBRK(L+1)........BREAK POINT SEQUENCE OF PP-REP. 
//    PCOEF(K,IPC,NCD).PP-COEFFICIENTS , I.E. 
//             PCOEF(J,I,,)=D**(J-1)(F(RBRK(I))   1<=I<=L. 
//    L................INDICATES NUM OF INTERVAL OF PP-REP. 
//    K................ORDER OF PP-REP(B-REP) 
//    NCD..............SPACE DIMENSION OF THE PP-REP 
//    N,T(N+K).........KNOT VECTOR OF THE BREP. N IS B-REP DIMENSION. 
//    IPC..............LENGTH OF 2ND DIMENSIONED ARRAY OF THE VARIALBLE 
//                     PCOEF. 
//    IPCW.............LENGTH OF 2ND DIMENSIONED ARRAY OF THE VARIALBLE 
//                     PCWORK. IPCW MUST BE GRETAER OR EQUAL TO (N-K+2). 
//    IRC..............ROW DIMENSION OF THE VARIABLE RCOEF 
// ***OUTPUT*** 
//    RCOEF(IRC,NCD)...B-COEFFICIENTS 
// ***WORK*** 
//    PCWORK(K,IPCW,NCD)........IS WORK ARRAY. IPCW>=(N-K+2). 
void blctpb_(
	const double *rbrk, const double *pcoef, int l,int k,
	int ncd, int n, const double *t, int ipc, int ipcw, int irc,
	double *pcwork, double *rcoef
);
