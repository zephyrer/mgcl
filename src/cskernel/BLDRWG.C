/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/bler.h"
#include "cskernel/bpval2.h"
#include "cskernel/bk1fli.h"
#include "cskernel/blcbp.h"
#include "cskernel/bldrw1.h"
#include "cskernel/bldrwg.h"

//     DIMENSION OF WK1 IS WK1(4*K*K+3*K) AS SPECIFIED IN COMMENT 
//   , BUT THE DECLARATION IS WK1(K,3) TO AVOID COMPILE ERROR 
//     AND TO UTILIZE THE AREA. 

//      SUBROUTINE TO DRAW LINE 
//    ***** POLYLINE DRAWING (GPL) VERSION ***** 
// *** INPUT  * 
//     GPL....IS A SUBROUTINE NAME TO DRAW POLYLINE. THE eING 
//            SEQUENCE IS AS: 
//              CALL GPL(N,X,Y), WHERE N:NUM OF POINT, X(N) AND 
//               Y(N) ARE THEIR COORDINATES, (X(I) Y(I)) FOR 1<=I<=N. 
//     NPY.....NUMBER OF POINT TO BE DRAWN FOR THE LENGTH WIND(4) 
//             ( FOR THE WINDOW LENGTH OF Y-COORDINATE ) . 
//             MINIMUM NPY =0, IN THIS CASE ONLY ONE SPAN BETWEEN KNOTS. 

//     WIND(4)..WINDOW SIZE TO CLIP: 
//           ( WIND(1),(2) ) IS THE CENTER (X,Y) OF THE WINDOW, AND 
//           ( WIND(3),(4) ) IS ( WIDTH,HEIGHT ) OF THE WINDOW. 
//           WIND(3)<=0. INDICATES CLIPPING IS NOT NECESSARY. 
//         EVEN WHEN WIND(3)<=0. WIND(4) IS NECESSARY TO INPUT FOR NPY. 
//     KLINI..SPECIFIES HOW INPUT DATA CORRESPONDS TO SCREEN COORDINATE; 
//           = 1 : (T,RCOEFX)      IS (X,Y) 
//           = 2 : (RCOEFX,T)      IS (X,Y) 
//           OTHERWISE : (RCOEFX,RCOEFY) IS (X,Y) OF THE SCREEN. 
//            WHEN KLIN=1,2 RCOEFY IS DUMMY ARGUMENT, NOT USED. 
//     K,N,T(N+K),RCOEFX(N),RCOEFY(N)......ARE B-REP TO DRAW: 
//            ORDER, B-REP DIMENSION, KNOT VECTOR, AND B-COEFFICIENTS. 
//     NWK2.....SPECIFIES LENGTH OF WK2 AS WK2(NWK2,2) 
//              NWK2 MUST BE .GE. N/2 AND RECOMMENDED .GE.N , 
//           USED TO PUT POSITIONAL DATA TO SUBROUTINE GPL. 
// *** WORK  * 
//     WK1(4*K*K+3*K),WK2(NWK2,2),RW(N) 
void bldrwg_(S_fp gpl, int npy, const double *wind, 
	int klini, int k, int n, const double *t,const double *rcoefx,
	const double *rcoefy, int nwk2, double *wk1, double *wk2, double *rw){
    // Local variables 
    int lend;
    int clip;
    double flmj;
    int klin;
    double vlen, wbrk[2], tnow;
    int jmkp1, i, j, l;
    double x, y, rleny, t2;
	int mj;
    double tc, dt, xm[2], ym[2];
    int npoint, ip1, is1, is2, np1;
    double flm;
    int lpp,nrw,nrwm1;
	const double* rcoef[2];

    if(nwk2<=2)
		return;
	rcoef[0]=rcoefx, rcoef[1]=rcoefy;

    // Parameter adjustments 
    --wind;
    wk1 -= k+1;
    --rw;
    --rcoefy;
    --rcoefx;
    --t;
    wk2 -= nwk2+1;

	// ***** INITIAL SET. 
    klin=klini;
    if(klin!=1 && klin!=2)
		klin = 0;
    
    rleny = (double) npy/wind[4];
    np1 = n+1;
    clip = wind[3]>0.f;

    if(clip){
		x = wind[3] * .5f;
		xm[0] = wind[1] - x;
		xm[1] = wind[1] + x;
		y = wind[4] * .5f;
		ym[0] = wind[2] - y;
		ym[1] = wind[2] + y;
// ***** OBTAIN INTERSECTION PARAM VAL WITH FRAME. 
		if (klin == 1) {
			bldrw1_(ym,xm,klin,k,n,&t[1],rcoef,
				&wk1[k+1],&wk2[nwk2+1],&nrw,&rw[1]);
		} else {
			bldrw1_(xm,ym,klin,k,n,&t[1],rcoef,
				&wk1[k+1],&wk2[nwk2+1],&nrw,&rw[1]);
		}
		if(nrw<=1)
			return;
	}else{
		nrw = 2;
		rw[1] = t[k];
		rw[2] = t[np1];
    }

// ****** CONVERT TO PP AND DRAW EACH LINE 
    nrwm1 = nrw-1;
    for(i=1; i<=nrwm1; ++i){
		ip1 = i+1;
		if(clip){
			tc = (rw[i]+rw[ip1])*.5f;
//   *** TEST IF CURRENT LINE IN FRAME OR NOT 
		    y = tc;
		    x = bler_(k, n, &t[1], &rcoefx[1], tc, 0);
			if (klin == 1) {
				y = x;
				x = tc;
		    }else if(klin == 0){
				y = bler_(k, n, &t[1], &rcoefy[1], tc, 0);
			}
		    if(x<=xm[0] || x>=xm[1] || y<=ym[0] || y>=ym[1])
				continue;
		}
// ***** DRAW LINE FROM RW(I) TO RW(I+1) 
		is1=bk1fli_(np1, &t[1], rw[i]);
		is2=bk1fli_(np1, &t[1], rw[ip1]);
		flm = rleny * (rw[ip1] - rw[i]) / (double) (is2 - is1 + 1);
//   *** COMPUTE FIRST POSITION'S DATA 
		t2 = rw[i];
		y = t2;
		x = bler_(k, n, &t[1], &rcoefx[1], t2, 0);
		if(klin == 1) {
			y = x;
		    x = t2;
		}else if(klin==0){
			y=bler_(k,n,&t[1],&rcoefy[1],t2,0);
		}
		npoint = 1;
		wk2[npoint + nwk2] = x;
		wk2[npoint + (nwk2 << 1)] = y;
//   *** DRAW LINE ONE SPAN BY ONE 
		for (j = is1; j <= is2; ++j) {
		    if(t[j]==t[j+1])
				continue;
//     ...CONVERT ONE SPAN TO PP-REP 
		    jmkp1 = j - k + 1;
			blcbp_(k,k,&t[jmkp1],&rcoefx[jmkp1],&wk1[k*3+1],wbrk,&wk1[k+1],&lpp);
		    vlen = fabs(wk1[k+2]);
		    if(klin==0){
				blcbp_(k,k,&t[jmkp1],&rcoefy[jmkp1],&wk1[k*3+1],wbrk,&wk1[(k<<1)+1],&lpp);
				vlen = sqrt(vlen*vlen+wk1[(k<<1)+2]*wk1[(k<<1)+2]);
			}
//     ...NOW PP-REP IN WK1(.,1-2) 
			tnow = t2;
		    t2 = min(rw[ip1],wbrk[1]);
		    mj = (int)(flm*vlen) + 2;
		    flmj = (double) mj;
			dt = (t2 - tnow) / flmj;
		    if (j == is2) {
			lend = mj - 1;
			} else {
			lend = mj;
		    }
//     ...DRAW LINE TO T2 BY INCREMENTING DT 
		    for (l = 1; l <= lend; ++l) {
			tnow += dt;
			y = tnow;
			x = bpval2_(wbrk, &wk1[k + 1], k, tnow, 1, 0);
			if (klin == 1) {
			    y = x;
				x = tnow;
			} else if (klin == 0) {
			    y = bpval2_(wbrk, &wk1[(k<<1)+1], k, tnow, 1, 0);
			}
			++npoint;
			wk2[npoint+nwk2] = x;
			wk2[npoint+(nwk2<<1)] = y;
			if(npoint>=nwk2){
			    (*gpl)(&npoint,&wk2[nwk2+1],&wk2[(nwk2<<1)+1]);
				wk2[nwk2+1] = wk2[npoint+nwk2];
			    wk2[(nwk2<<1)+1] = wk2[npoint+(nwk2<<1)];
			    npoint=1;
			}
		    }
		}
// DRAW LAST 1-SPAN. 
		tnow = t2;
		y = tnow;
		x = bpval2_(wbrk, &wk1[k+1], k, tnow, 1, 0);
		if (klin == 1) {
			y = x;
		    x = tnow;
		} else if (klin == 0) {
			y = bpval2_(wbrk, &wk1[(k<<1)+1], k, tnow, 1, 0);
		}
		++npoint;
		wk2[npoint+nwk2]=x;
		wk2[npoint+(nwk2<<1)]=y;
		(*gpl)(&npoint, &wk2[nwk2 + 1], &wk2[(nwk2 << 1) + 1]);
    }
}
