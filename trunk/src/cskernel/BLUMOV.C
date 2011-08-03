/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <malloc.h>
#include "cskernel/tolerance.h"
#include "cskernel/Ble.h"
#include "cskernel/blumor.h"
#include "cskernel/b2vb.h"
#include "cskernel/bk2fli.h"

// BLUMOV MOVES A LINE B-REP AS SPECIFIED IN KT. 
// *** INPUT * 
//       K,N,T(N+K),RCOEFI(IRCI,NCD),IRCI,NCD....DESCRIBE B-REP INPUT 
//             K:ORDER,     NCD:SPACE DIMENSIONN,   N:B-REP DIMENSION 
//                T(N+K):KNOT VECTOR,  RCOEFI(.,.): B-COEFFICIENTS 
//       TC...IS THE CENTER OF TRANSLATION, PARAMETER VALUE OF THE B-REP
//       P(NCD).....IS THE NEW POINT COORDINATE, 
//                THE POINT F(TC) IS MOVED TO P(.), WHERE F IS B-REP 
//       KT.....INDICATES WHAT KIND OF MOVE IS BEING PERFORMED. 
//          =1 : TWO END POINTS FIXED AND CENTER OF MOVE BETWEEN THEM. 
//          =2 : ONE POINT,TFI(1),FIXED. THE OTHER POINT IS FREE 
//          =3 : TWO POINTS,TFI(.),FIXED AND CENTER OF MOVE BETWN THEM. 
//          =4 : NO FIXED POINT SPECIFIED, MINIMUM MOVE 
//          =5 : PARALLEL TRANSLATION OF THE LINE 
//       TFI(2)....SPECIFIES PARAMETER VALUES OF B-REP TO FIX AT. 
//          WHEN KT=2, TFI(2) IS DUMMY. 
//          WHEN KT=1, 4 AND 5, TFI(.) ARE DUMMY. 
//       IRC.....ROW DIMENSIONS OF THE VARIABLES RCOEF 
// *** OUTPUT * 
//       RCOEF(IRC,NCD).....ARE NEW TRANSLATED B-COEF'S. 
//          RCOEF MAY BE THE SAME AREA AS RCOEFI 
//       TFO(2).....ARE PARAMETER VALUES OF THE B-REP, GIVES TWO 
//          BOUNDARY POINTS OF THE MOVE. I.E., MOVE IS PERFORMED IN 
//          TFO(1) < T < TFO(2). 
// ***WORK* 
//       WRATIO(N)...WORK AREA OF LENGTH N. 
void blumov_(int k, int n, const double *t, 
	const double *rcoefi, int irci, int ncd, double tc, 
	const double *p, int kt,const double *tfi, int irc, double *wratio,
	double *rcoef, double *tfo)
{
    int ism1, istrt, iend, i, j, l, m, jj, kmov, np1, lmk, idx;
    double xsum, ratio, te, ts, tk, tnp1;
	double delx_buf[20];
    double x_buf[20];
	double* delx=delx_buf;
	double* x=x_buf;

    // Parameter adjustments 
    --wratio;
    --p;
    rcoefi -= irci+1;
    rcoef -= irc+1;

    np1 = n+1;
	tk=t[k-1], tnp1=t[n];//Start and end parameter values.
    if(tc<tk || tc>tnp1){
		for(m=1; m<=n; ++m){
			for(jj=1; jj<=ncd; ++jj)
				rcoef[m+jj*irc] = rcoefi[m+jj*irci];
	    }
	    tfo[0] = tk;
		tfo[1] = tnp1;
		return;
    }

	if(k>20){
		delx=(double*)malloc(sizeof(double)*(k));
		x=(double*)malloc(sizeof(double)*(k));
	}

	ble_(k,n,t,&rcoefi[irci+1],irci,ncd,tc,0,delx);
    for(jj=1; jj<=ncd; ++jj)
		delx[jj-1] = p[jj]-delx[jj-1];
    l=bk2fli_(np1,t,tc);
    lmk = l-k;
    b2vb_(k,n+k,t,tc,l,x);
//     SET I AND J FOR THE CASE OF KT=4 
    kmov = 1;
    i = l;
    j = l;

	ts = tfi[0];
	te = tfi[1];
	if(kt==5){
	// ====CASE OF NO FIXED POINTS, I.E. PARALLEL TRANSLATION. ==== 
	    for(m=1; m<=n; ++m){
			for(jj=1; jj<=ncd; ++jj)
				rcoef[m+jj*irc] = rcoefi[m+jj*irci]+delx[jj-1];
	    }
		istrt = 1;
	    iend = n;
	}else{
		if(kt==1){
		    i = k;
			j = n;
		}else if(kt==2){
		// ====CASE OF KT=2, ONE POINT FIXED.==== 
			if(tfi[0]<=tc){
				kmov=2;
				if(ts<tk)
				    ts=tk;
				j=n;
			}else{
				kmov = 3;
				te = tfi[0];
				if(te>tnp1)
				    te = tnp1;
				i=1;
			}
		}else if(kt==3){
		// ====CASE OF KT=3, TWO SPECIFIED POINTS ARE FIXED.==== 
    		if(ts>te){
				ts = te;
				te = tfi[0];
    		}
    		if(ts<tk)
				ts = tk;
    		if(te>tnp1)
				te = tnp1;
    		if(tc<ts){
				kmov = 3;
				te = ts;
				i = 1;
    		}else if(tc>te){
				kmov = 2;
				ts = te;
				j = n;
			}		
		}
		if(kt==2 || kt==3){
		// CASE OF KT=2,3 AND KMOV=1,2,3
		    if(kmov==1 || kmov==2)
				i=bk2fli_(np1,t,ts);
		    if(kmov==1 || kmov==3)
				j=bk2fli_(np1,t,te);
		}

		// GET TRANSLATION RATIO IN WRATIO(.) 
		blumor_(kmov,&i,&j,tc,k,n,t,&istrt,&iend,&wratio[1]);
		xsum = 0.f;
		idx = 1;
		for(m=lmk+1; m<=l; ++m){
			if(m>=istrt && m<=iend)
				xsum += x[idx-1]*wratio[m];
			++idx;
		}
		if(xsum<=bzmzro_()){
			for(m=1; m<=n; ++m){
				for(jj=1; jj<=ncd; ++jj)
					rcoef[m+jj*irc] = rcoefi[m+jj*irci];
			}
			goto L7777;
		}
		//    GET NEW B-COEF'S 
		ism1 = istrt-1;
		for(m=1; m<=ism1; ++m){
			for(jj=1; jj<=ncd; ++jj){
				rcoef[m+jj*irc] = rcoefi[m+jj*irci];
			}
		}
		for(m=istrt; m<=iend; ++m){
			ratio = wratio[m]/xsum;
			for(jj=1; jj<=ncd; ++jj)
				rcoef[m+jj*irc] = rcoefi[m+jj*irci]+delx[jj-1]*ratio;
		}
		for(m=iend+1; m<=n; ++m){
			for(jj=1; jj<=ncd; ++jj)
				rcoef[m+jj*irc] = rcoefi[m+jj*irci];
		}
	}

    tfo[0] = tc;
    tfo[1] = tc;
    idx = istrt;
    if(idx<k)
		idx = k;
    tfo[0] = t[idx-1];
    idx = iend+k;
    if(idx>np1)
		idx = np1;
    tfo[1] = t[idx-1];
L7777:
	if(k>20) {free(delx); free(x);}
}
