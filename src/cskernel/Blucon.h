/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//    BLUCON WILL CONNECT TWO B-REP. CURVES AND GET ONE B-REP. 
// *** INPUT  * 
//     K,N1,T1(N1+K),RCOEF1(IRC1,NCD),IRC1,NCD....1ST B-REP 
//      *** N1 MAY BE 0, IN THIS CASE B-REP 2 COPIED INTO (N1,T1,RCOEF1) 
//     N2,T2(N2+K),RCOEF2(IRC2,NCD),IRC2...2ND B-REP TO CONCATENATE. 
//     JCNI.....C(JCNI) CLASS ABOUT CONTINUITY 
// *** OUTPUT * 
//     N1,T1(N1+K),RCOEF1(.,.).....NEW B-REP CONCATENATED. 
//     RATIO....RATIO OF PARAMETER VALUE REGION (NEW/OLD) 
//                                 OF 2ND B-REP. IN NEW B-REP. 
//     IT2S.....STARTING POINT OF 2ND B-REP. IN NEW B-REP. 
// *** WORK   * 
//     WORK(K,K),PCOEF(K,KK)....
//				WORK ARRAY OF EACH LENGTH, WHERE KK=MAX(NCD,K+2). 
void blucon_(int k, int *n1, double *t1, double *rcoef1, int irc1,
		int ncd, int n2, const double *t2, const double *rcoef2, int irc2, int jcni, 
		double *work, double *pcoef, double *ratio, int *it2s);
