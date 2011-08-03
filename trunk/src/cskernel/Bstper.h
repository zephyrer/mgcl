/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BSTPER COMPUTES TANGENT PLANE OF ONE OF GIVEN SURFACE B-REP PERIMETER. 
// *** INPUT *** 
// KU,LUD,UKT(LUD+KU),KV,LVD,VKT(LVD+KV),SURF(ISR1,ISR2,3)...SURFACE B-REP .
//  KR......PERIMETER NUM:    1     2     3     4 
//                          V=MIN  U=MAX V=MAX U=MIN 
// *** OUTPUT *** 
//  K,N,T(N+K),RCOEF(IRC,3).....T.P. B-REP OBTAINED (ORDER K=MAX(KU,KV)) 
// *** WORK ARRAY *** 
//  WORK(N,K*2): WORK AREA 
void bstper_(int ku, int lud,const double *ukt, 
	int kv, int lvd,const double *vkt,const double *surf, int isr1, int isr2, int ncd,
	int kr, int irc, double *work, int *k, int *n, double *t, double *rcoef);
