/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/b1nk.h"
#include "cskernel/bk2fli.h"

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
	double te, double *tw, double *batj, double *rbound)
{
    // System generated locals 
    int n2mk;

    // Local variables 
    int iemk, ismk, i, j, k, l, n2mkp1;
    double rcnew;
    int i2, n2;
    int ie, is, kp1, np1;

// **************************** START OF BLQBOX  ************************

    // Parameter adjustments 
	k=kp;
    kp1 = k + 1;
    batj -= kp1;
    --t;
    rbound -= 3;
    rcoef -= irc + 1;;
    --tw;

    // Function Body 
    np1 = n + 1;
    is=bk2fli_(np1, &t[1], ts);
    ie=bk2fli_(np1, &t[1], te);
    ismk = is - k;
    n2 = ie - ismk;
    n2mkp1 = n2 - k + 1;

//  <<<<< 1. Build new knot configuration. >>>>> 
    for(i=1; i<=k ; ++i) tw[i] = ts;
    if(is<ie)
		for(i=kp1; i<=n2; ++i) tw[i] = t[ismk+i];
    for(i=1; i<=k; ++i) tw[n2+i] = te;
//  <<<<< 2. Compute Minimum and Maximun of NEW B-COEFFICIENTS.>>>>> 
//   .....2.1 Starting K Coefficients. 
    b1nk_(kp, &t[1], is, &tw[1], 1, &batj[kp1]);
    for(j=1; j<=ncd; ++j){
		rcnew = 0.f;
		for(l=1; l<=k; ++l)
		    rcnew += batj[l+k*k]*rcoef[ismk+l+j*irc];
		rbound[(j<<1) + 1] = rcnew;
		rbound[(j<<1) + 2] = rcnew;
		}
    for(i=2; i<=k; ++i){
		b1nk_(kp, &t[1], is, &tw[1], i, &batj[kp1]);
		for (j = 1; j <= ncd; ++j) {
			rcnew = 0.f;
		    for(l=1; l<=k; ++l)
				rcnew += batj[l+k*k]*rcoef[ismk+l+j*irc];
			if(rbound[(j<<1) + 1] > rcnew) rbound[(j<<1) + 1] = rcnew;
			if(rbound[(j<<1) + 2] < rcnew) rbound[(j<<1) + 2] = rcnew;
		}
    }
//   .....2.2 Coefficients  After Starting K Coefficients 
//            and Before Ending K Coefficients 
    if (n2mkp1 > kp1) {
	n2mk = n2 - k;
	for(i=kp1; i<=n2mk; ++i){
		for(j=1; j<=ncd; ++j){
			if(rbound[(j<<1)+1] > rcoef[ismk+i+j*irc])
				rbound[(j<<1)+1] = rcoef[ismk+i+j*irc];
			if(rbound[(j<<1)+2] < rcoef[ismk+i+j*irc])
				rbound[(j<<1)+2] = rcoef[ismk+i+j*irc];
		}
	}
    }
//   .....2.3 Ending K Coefficients. 
    i2 = n2mkp1;
    if(i2<kp1) i2 = kp1;
    if(i2<=n2){
	iemk = ie-k;
	for(i=i2; i<=n2; ++i){
	    b1nk_(kp, &t[1], ie, &tw[1], i, &batj[kp1]);
	    for(j=1; j<=ncd; ++j){
		rcnew = 0.f;
		for(l=1; l<=k; ++l)
			rcnew += batj[l+k*k]*rcoef[iemk+l+j*irc];
		if(rbound[(j<<1)+1] > rcnew) rbound[(j<<1)+1] = rcnew;
		if(rbound[(j<<1)+2] < rcnew) rbound[(j<<1)+2] = rcnew;
	    }
	}
    }
}
