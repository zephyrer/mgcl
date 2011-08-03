/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/bkcrng.h"

//  BLUREV REVERSES THE DIRECTION OF PARAMETER OF A GIVEN B-REP. 
// *** INPUT  * 
//      K,N1,T1(N1+K),RCOEF1(IRC1,NCD),IRC1,NCD....ORIGINAL B-REP. 
//      IRC2........ROW DIMENSION OF RCOEF2 
// *** OUTPUT * 
//      N2,T2(N2+K),RCOEF2(IRC2,NCD)...NEW B-REP OBTAINED. N2=N1 
// *** NOTE * 
//      N2,T2,RCOEF2 MAY BE THE SAME AREA AS N1,T1,RCOEF1, EACH. 
void blurev_(int k, int n1, const double *t1, 
	const double *rcoef1, int irc1, int ncd, int irc2, 
	int *n2, double *t2, double *rcoef2)
{
    // Local variables 
    int n1by2, j;
    double rsave;
    int i1, i2, npkby2;
    double te;
    double ts;
    // Parameter adjustments 
    --t2;
    --t1;
    rcoef1 -= irc1+1;
    rcoef2 -=  irc2+1;;

    *n2 = n1;
// Change the B-coefficients of the B-REP. 
    n1by2 = n1/2;
    for(j=1; j<=ncd; ++j){
		i2 = n1;
		for(i1=1; i1<=n1by2; ++i1){
			rsave = rcoef1[i2+j*irc1];
		    rcoef2[i2+j*irc2] = rcoef1[i1+j*irc1];
		    rcoef2[i1+j*irc2] = rsave;
			--i2;
		}
		if(i2>n1by2)
			rcoef2[i2+j*irc2] = rcoef1[i2+j*irc1];
    }
// *** CHANGE OF KNOT VECTOR *** 
//  1. Reverse the Knot configuration. 
    ts = t1[k];
    te = t1[n1+1];
    i2 = n1+k;
    npkby2 = i2/2;
    for(i1=1; i1<=npkby2; ++i1){
		rsave = t1[i1];
		t2[i1] = -t1[i2];
		t2[i2] = -rsave;
		--i2;
    }
    if(i2>npkby2)
		t2[i2] = -t1[i2];
//  2. Parameter Range adjustment. 
    bkcrng_(k,*n2,&t2[1],ts,te,&t2[1]);
}
