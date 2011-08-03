/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLCPB CONVERTS PP-REP(RBRK,PCOEF,L) TO B-REP(T,RCOEF,N). 
// THE PP-REP MUST HAVE THE CONTINUITY K-2 AT EACH BREAK POINT. 
// ***INPUT* 
//    RBRK(L+1)........BREAK POINT SEQUENCE OF PP-REP. 
//    PCOEF(K,IPC,NCD).PP-COEFFICIENTS , I.E. 
//             PCOEF(J,I,,)=D**(J-1)(F(RBRK(I))   1<=I<=L. 
//    L................INDICATES NUM OF INTERVAL OF PP-REP. 
//    K................ORDER OF PP-REP(B-REP) 
//    NCD..............SPACE DIMENSION OF THE PP-REP 
//    IPC..............LENGTH OF 2ND DIMENSIONED ARRAY OF THE VARIALBLE 
//                     PCOEF, MUST BE .GE.L+1 . 
//    IRC..............ROW DIMENSION OF THE VARIABLE RCOEF 
//    IT...............INDICATES WHETHER KNOT VECTOR T(.) IS INPUT 
//                     OR NOT(GENERATED IN THE PREVIOUS CALL OF BLCPB) 
//                        =1  BLCPB GENERATES T 
//                        =2  T IS INPUT(BLCPB NOT GENERATE T) 
//    T(N+K)...........KNOT VECTOR (NECESSARY WHEN IT=2) 
// ***OUTPUT*** 
//    T(N+K)...........KNOT VECTOR OF B-REP(WHEN IT=1) 
//    RCOEF(IRC,NCD)...B-COEFFICIENTS 
//    N................B-REP DIMENSION 
void blcpb_(
	const double *rbrk, double *pcoef, int l, int k, int ncd, int ipc, int irc, int it, 
	double *t, double *rcoef, int *n
);
