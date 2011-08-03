/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/blurev.h"
#include "cskernel/bkcrng.h"

// BLUCPR CHANGES PARAMETER RANGE OF GIVEN B-REP, MAY REVERSE THE PARAM 
// DIRECTION. 
// ***INPUT*** 
//   S1,S2.....GIVE NEW PARAM RANGE, S1 CORRESPONDS TO T1(K) AND 
//             S2 TO T1(N+1). SO IF S1 > S2, GIVEN B-REP'S PARAMETER 
//             DIRECTION SHOULD BE REVERSED. 
//   K,N1,T1(N+K),RCOEF1(IRC1,NCD),IRC1,NCD.....PROVIDE ORIGINAL B-REP 
//             ORDER,SPACE DIMENSION, B-REP DIMENSION, KNOT VECTOR, 
//             B-COEF'S, EACH. 
//   IRC2......ROW DIMENSION OF THE VARIABLE RCOEF2. 
// ***OUTPUT*** 
//   N2,T2(N2+K),RCOEF2(IRC2,NCD).....UPDATED B-REP 
// *** NOTE*** 
//   AREA (N1,T1,RCOEF1) AND (N2,T2,RCOEF2) MAY BE THE SAME, EACH. 
void blucpr_(double s1, double s2, int k, 
	int n1, const double *t1, const double *rcoef1, int irc1, 
	int ncd, int irc2, int *n2, double *t2, double *rcoef2)
{
    int i, j;
    double te, ts;

    if(s1>=s2){
		ts = s2;
		te = s1;
		blurev_(k,n1,t1,rcoef1,irc1,ncd,irc2,n2,t2,rcoef2);
	    bkcrng_(k, n1, t2, ts, te, t2);
		return;
    }

    // Parameter adjustments 
    rcoef1 -= irc1+1;
	rcoef2 -= irc2+1;

	ts = s1;
	te = s2;
	for(i=1; i<=n1; ++i){
	    for(j=1; j<=ncd; ++j)
			rcoef2[i+j*irc2] = rcoef1[i+j*irc1];
	}
	*n2 = n1;
    bkcrng_(k, n1, t1, ts, te, t2);
}
