/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/b1nk.h"
#include "cskernel/bk2fli.h"

// BLUNK GENERATES B-COEFFICIENTS OF NEW KNOT CONFIGURATION T2. 
// *** INPUT * 
//     K,N1,T1(N1+K),RCOE1(IRC1,NCD),IRC1,NCD...DESCRIBE THE OLD B-REP 
//     N2,T2(N2+K).....KNOT VECTOR OF THE NEW B-REP 
//     IRC2....ROW DIMENSION OF THE VARIABLE RCOE2 
// *** OUTPUT * 
//     RCOE2(IRC2,NCD)..THE NEW B-COEF OBTAINED 
// *** WORK * 
//     BATJ(K,K)  LENGTH OF K*K 
// *** NOTE * 
//     THE NEW KNOT T2 MUST SATISFY THE FOLLOWLING CONDITION; 
//          FOR ANY I (1<=I<=N2-1), THE NUM OF J S.T. 
//          T2(I)<= T1(J) < T2(I+1) 
//          MUST BE LESS OR EQUAL 1. 
void blunk_(int k, int n1, const double *t1, 
	const double *rcoe1, int irc1, int ncd, int n2, 
	const double *t2, int irc2, double *batj, double *rcoe2)
{
    // System generated locals 
    int batj_offset;

    // Local variables 
    int mumk, i, j, l, isame;
    double rcnew;
    int i1=0;
    int mu, km1, n1p1;

// **************************** START OF BLUNK  ************************
// GET NEW B-COEFFICIENTS. 
    // Parameter adjustments 
    batj_offset = k+1;
    batj -= batj_offset;
    --t1;
    rcoe1 -= irc1+1;
    --t2;
    rcoe2 -= irc2+1;

    // Function Body 
    km1 = k-1;
    n1p1 = n1+1;
    isame = 0;
//   ISAME IS THE FLAG TO INDICATE THE PREVIOUS COEF IS OBRAINED BY 
//   NO COMPUTATION, SINCE THE NEW KNOT THE SAME KNOT CONFIGURATION. 
//            ISAME=0: NOT THE SAME, =1:THE SAME KNOTS. 
    for(i=1; i<=n2; ++i){
		if(isame==1) {
			++i1;
			if(t2[i+km1] == t1[i1+km1]){
				for(j=1; j<=ncd; ++j)
					rcoe2[i+j*irc2] = rcoe1[i1+j*irc1];
				continue;
			}
		}
		mu=bk2fli_(n1p1, &t1[1], t2[i]);
		if(isame!=1){
		// ===== 1. TEST IF THE SAME KNOT COFIGURATION ===== 
			//    .... (1) TAKE THE KNOT-MULTIPLICITY INTO ACCOUNT. 
			for(l=1; l<=k; ++l)
				if(t2[i+l] != t2[i])
					break;
			i1 = mu-l+1;
			//     .... (2) TEST OF SAME KNOT COFIG. 
			for(l=1; l<=km1; ++l){
				if(t2[i+l] != t1[i1+l])
					goto L40;
			}
			//     .... (3) NOW THE SAME KNOTS, AND SO NO COMPUTATION 
			for(j=1; j<=ncd; ++j)
				rcoe2[i+j*irc2] = rcoe1[i1+j*irc1];
			isame=1;
			continue;
		}
L40:
		// ===== 2. THE CASE THAT COMPUTATION BY OTHLO ALGORITHM IS NECESSARY == 
		mumk = mu-k;
		b1nk_(k,&t1[1],mu,&t2[1],i,&batj[batj_offset]);
		for(j=1; j<=ncd; ++j){
			rcnew = 0.f;
		    for(l=1; l<=k; ++l){
				i1 =mumk+l;
				if(i1<1)
					i1 = 1;
				if(i1>n1)
					i1 = n1;
				rcnew += batj[l+k*k]*rcoe1[i1+j*irc1]
				;
		    }
		    rcoe2[i+j*irc2] = rcnew;
		}
		isame = 0;
    }
}
