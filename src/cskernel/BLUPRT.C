/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/Blunk.h"
#include "cskernel/bk2fli.h"

// BLUPRT computes a partial b-spline, given a spline and 
// a partial parameter range of the original knot configuration. 
// The new spline is exactly the same as the old one, 
// although the representation is partial. 
// *** INPUT * 
//     K,N1,T1(N1+K),RCOE1(IRC1,NCD),IRC1,NCD... 
//                     Describe the original B-REP. 
//     TSI,TEI.....New parameter range, should be 
//               TSI < TEI.  If TSI <= T1(K) or T1(N1+1) <= TEI, 
//               obtained b-rep is not partial and the same as 
//               the original one, regarding to the start or end 
//               side of the origial line each. 
//     IRC2....ROW DIMENSION OF THE VARIABLE RCOE2 
//     multiple...indicates if knot multiplicity of K is necessary at start
//                and end parameter of T2.
//               =0: unnecessary, !=0: necessary.
// *** OUTPUT * 
//     N2,T2(N2+K),RCOE2(IRC2,NCD)..THE NEW B-COEF OBTAINED. 
// *** WORK * 
//     WORK(K,K)  LENGTH OF K*K 
// *** NOTE * 
//     If TSI >= TEI, N2=0 is set and nothing will be done. 
void bluprt_(int k, int n1, const double *t1, 
	const double *rcoe1, int irc1, int ncd, double tsi, 
	double tei, int irc2, double *work, int *n2, 
	double *t2, double *rcoe2, int multiple)
{
    // System generated locals 
    int rcoe1_offset, 
	    rcoe2_offset;

    // Local variables 
    int ismk, i, j;
    double ts,te;
    int is,ie,n1p1;

    // Parameter adjustments 
    --t1;
    rcoe1_offset = irc1+1;
    rcoe1 -= rcoe1_offset;
    rcoe2_offset = irc2+1;
    rcoe2 -= rcoe2_offset;
    --t2;

    // Function Body 
    *n2 = 0;
    if(tsi>=tei)
		return;

    n1p1 = n1+1;
    ts=t1[k];
    if(ts<tsi)
		ts = tsi;
    te=t1[n1p1];
    if(te>tei)
		te = tei;
    if(ts>=te)
		return;

    is=bk2fli_(n1p1,&t1[1],ts);
    ie=bk2fli_(n1p1,&t1[1],te);
    while(te == t1[ie])
		--ie;
    ismk = is-k;
    *n2 = ie-ismk;

//  Setup the new kont confuguration. 
// ... 1. Starting knot configuration 
    if(!multiple && ts==t1[is]){
		//Case that starting parameter is the same as one of original knots.
		for(i=1; i<=k; ++i)
			t2[i] = t1[ismk+i];
	}else{
		for(i=1; i<=k; ++i)
			t2[i] = ts;
    }
// ... 2. Midpoint knot confuguration. 
    for(i=k+1; i<=*n2; ++i)
		t2[i] = t1[ismk+i];
// ... 2. Ending point knot confuguration. 
    if(!multiple && te==t1[ie+1]){
		//Case that end point IS the same as one of the knots.
		for(i=1; i<=k; ++i)
			t2[*n2+i] = t1[ie+i];
    }else{
		for(i=1; i<=k; ++i)
			t2[*n2+i] = te;
    }

// Get new B-Coefficients. 
    if(ts==t1[is] && te==t1[ie+1]){
		for(i=1; i<=*n2; ++i){
		    for(j=1; j<=ncd; ++j)
				rcoe2[i+j*irc2] = rcoe1[ismk+i+j*irc1];
		}
    }else{
		blunk_(k,n1,&t1[1],&rcoe1[rcoe1_offset],irc1,ncd,*n2,&t2[1], 
			irc2,work,&rcoe2[rcoe2_offset]);
    }
}
