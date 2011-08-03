/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// BLUAKT ADDS KNOT(S) INTO ORIGINAL B-REP, GUARANTREEING THE SAME LINE.
// *** INPUT  * 
//     K,N1,T1(N1+K),RCOEF1(IRC1,NCD),IRC1,NCD..... ORIGINAL B-REP. 
//     NAD,TAD(NAD),MLT(NAD).....NUM OF KNOT, PARAMETER VALUE, AND 
//              MULTIPLICITY AT THE PARAM TAD(.). 
//     IRC2....ROW DIMENSION OF RCOEF2 
// *** OUTPUT * 
//     N2,T2(N2+K),RCOEF2(IRC2,NCD)....NEW B-REP OBTAINED. 
// *** WORK   *   WORK(K,K) 
void bluakt_(int k, int n1, const double *t1, 
	const double *rcoef1, int irc1, int ncd, int nad, 
	const double *tad, const int *mlt, int irc2, double *work, 
	int *n2, double *t2, double *rcoef2
);
