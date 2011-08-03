/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include "cskernel/tolerance.h"
#include "cskernel/blipp.h"
#include "cskernel/bk1fli.h"
#include "cskernel/bler.h"

// BLDRW2 IS LOCAL SUB FOR BLDRW1, STORES PARAM VALUE TT SO THAT 
// INCREASING PROPERTY OF RW(.) HOLDS. 
void bldrw2_(int *nrw, double *rw, double tt, double terror){
    int is,i;

    if(*nrw >= 1){
		is=bk1fli_(*nrw, rw, tt);
		if(is>=1){
			if(tt-terror <= rw[is-1])
				return;
		}
		++is;
		if(is<=*nrw) {
			if(tt+terror>=rw[is-1])
				return;
		}
		for(i=*nrw; i>=is; --i)
			rw[i]=rw[i-1];
		rw[is-1] = tt;
		++(*nrw);
    }else{
		*nrw = 1;
		rw[0] = tt;
    }
}

// INTERNAL SUBROUTINE FOR BLDRWC AND BLDRWG, GETS INTERSECTION PARAM 
// VALUES WITH GIVEN FRAME (XM(2),YM(2)) 
// *** INPUT  * 
//     XM(2),YM(2)....WINDOW TO CLIP: 
//         ( XM(1),(2) ) ARE MIN AND MAX FOR RCOEF[0] 
//         ( YM(1),(2) ) ARE MIN AND MAX FOR RCOEF[1]  (KLIN.EQ.0) 
//                                       FOR KNOT T (KLIN.NE.0) 
//     KLIN...SPECIFIES HOW INPUT DATA CORRESPONDS TO FRAME (XM,YM); 
//           .NE.0 : (T,RCOEF[0])     IS (X,Y) 
//           .EQ.0 : (RCOEF[0],RCOEF[1]) IS (X,Y) 
//     K,N,T(N+K),RCOEF[0](N),RCOEF[1](N)......ARE B-REP TO DRAW. 
// *** OUTPUT * 
//     NRW,RW(NRW)......ARE INTERSECTION PARAM VALUES WITH FRAME. 
// *** WORK *       W1(4*K*K+3*K),W2(N) 
void bldrw1_(const double *xm,const double *ym, int klin, 
	int k, int n, const double *t, const double **rcoef, 
	double *w1, double *w2, int *nrw, double *rw)
{
    // Local variables 
    double terror, xerror, yerror;
    int i, l, iend;
    double x, y;
    double te, rzero;
    int nx, ny;
    double ts, xn[2], yn[2], tsp;

// ***** START OF EXECUTION 
    // Parameter adjustments 
    --xm;
    --ym;
    --w2;
    --t;

    // Function Body 
    ts = t[k];
    te = t[n+1];
    tsp = ts-1.f;
// Computing MAX 
	rzero=bzrzro_();
    terror = (te - ts) * rzero;
    xerror = (xm[2] - xm[1]) * rzero;
    yerror = (ym[2] - ym[1]) * rzero;
    xn[0] = xm[1] + xerror;
    xn[1] = xm[2] - xerror;
    yn[0] = ym[1] + yerror;
    yn[1] = ym[2] - yerror;
    *nrw = 0;
// ***** START POINT CHECK 
    x = bler_(k,n,&t[1],rcoef[0],ts,0);
    y = ts;
    if(klin == 0)
		y=bler_(k, n, &t[1], rcoef[1], ts, 0);
    if(x >= xn[0] && x <= xn[1] && (y >= yn[0] && y <= yn[1])) {
		*nrw = 1;
		rw[0] = ts;
		tsp = ts;
    }
// ***** OBTAIN INTRSCT POINT WITH X=X-MIN AND X-MAX 
    for(l=1; l<=2; ++l){
		blipp_(k,n,&t[1],rcoef[0],xm[l],xerror,tsp,n,w1,&nx,&w2[1],&iend);
		for(i=1; i<=nx; ++i){
		    y = w2[i];
			if(klin == 0)
				y=bler_(k, n, &t[1], rcoef[1], y, 0);
		    if(y >= yn[0] && y <= yn[1])
				bldrw2_(nrw, rw, w2[i], terror);
		}
    }
// ***** OBTAIN INTRSCT POINT WITH Y=Y-MIN AND Y-MAX 
    for(l=1; l<=2; ++l){
		if(klin==0){
		    blipp_(k,n,&t[1],rcoef[1],ym[l],yerror,tsp,n,w1,&ny,&w2[1],&iend);
		    for(i=1; i<=ny; ++i){
				x=bler_(k,n,&t[1],rcoef[0],w2[i],0);
				if(x>=xn[0] && x<=xn[1])
					bldrw2_(nrw, rw, w2[i], terror);
			}
		}else if(ym[l]>=ts && ym[l]<=te){
			x=bler_(k,n,&t[1],rcoef[0],ym[l],0);
			if(x>=xn[0] && x<=xn[1])
				bldrw2_(nrw, rw, ym[l], terror);
		}
    }
// ***** END POINT CHECK 
    x=bler_(k,n,&t[1],rcoef[0],te,0);
    y=te;
    if(klin == 0)
		y=bler_(k, n, &t[1], rcoef[1], te, 0);
    if(x>=xn[0] && x<=xn[1] && (y>=yn[0] && y<=yn[1]))
		bldrw2_(nrw, rw, te, terror);
}
