/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
#include <math.h>
#include "cskernel/tolerance.h"
#include "cskernel/bler.h"
#include "cskernel/bpval2.h"
#include "cskernel/bk2fli.h"
#include "cskernel/blcbp.h"
#include "cskernel/bldrwcr.h"

///****Rational version of bldrwc_.*****
//      SUBROUTINE TO DRAW LINE AND LINE NAME 
// *** INPUT  * 
//     kfunc .......is function kind of MOVEA2,LINEA2:
//          1:         movea2(int x, int y);
//			2:         movea2(float x, float y);
//			otherwise: movea2(double x, double y);
//			regarding to linea2, the same.							
//     MOVEA2,LINEA2....ARE SUBROUTINE NAMES TO MOVE CURRENT POSITION 
//           AND, TO DRAW (STRAIGHT LINE-SEGMENT) FROM CURRENT POSITION 
//           TO SPECIFIED POSITION; THEIR CALLING SEQUENCE ARE AS: 
// CALL MOVEA2(X,Y)....MOVE PEN POSITION TO (X,Y) WITHOUT DRAWING, 
// CALL LINEA2(X,Y)....DRAW LINE-SEGMENT FROM CURRENT POSITION TO (X,Y). 
//     NPY.....NUMBER OF POINT TO BE DRAWN FOR THE LENGTH WIND(4) 
//             ( FOR THE WINDOW LENGTH OF Y-COORDINATE ) . 
//     WIND(4)..WINDOW SIZE TO CLIP: 
//           ( WIND(1),(2) ) IS THE CENTER (X,Y) OF THE WINDOW, AND 
//           ( WIND(3),(4) ) IS ( WIDTH,HEIGHT ) OF THE WINDOW. 
//           WIND(3)<=0. INDICATES CLIPPING IS NOT NECESSARY. 
//         EVEN WHEN WIND(3)<=0. WIND(4) IS NECESSARY TO INPUT FOR NPY. 
//		nrw, rw[nrw]... are parameter ranges of after clipping.
//          generally rw[i] are intersection point parameters with x=minimum,
//			x=maximum , y=minimum, and y=maximum of the clipping window.
//			rw[i] must be increading order for 0<=i<=nrw-1.
//			*** When WIND(3)<=0 (clipping is unnecessary), nrw and rw are not
//			input, and used as work array for rw[0] and rw[1].
//			That is, these values will be destroyed.	                 
//     KLINI..SPECIFIES HOW INPUT DATA CORRESPONDS TO SCREEN COORDINATE; 
//           = 1 : (T,RCOEF[0])      IS (X,Y) 
//           = 2 : (RCOEF[0],T)      IS (X,Y) 
//                FOR KLINI=1,2, RCOEF[1] IS THE WEIGHT COEFFICIENTS. 
//           OTHERWISE : (RCOEF[0],RCOEF[1]) IS (X,Y) OF THE SCREEN. 
//                IN THIS CASE, RCOEF[2] IS THE WEIGHT COEFFICIENTS. 
//     K,N,T(N+K),rcoef[.][.]......ARE B-REP TO DRAW: 
//            ORDER, B-REP DIMENSION, KNOT VECTOR, AND B-COEFFICIENTS. 
// *** WORK  * 
//     WK1(3K+K*K) 
void bldrwcr_(int kfunc, S_fp movea2, S_fp linea2, unsigned npy, const double *wind,
			int nrw, double *rw, int klini, int k, int n,
			const double *t,const double **rcoef, double *wk1){
    int clip;
    int klin, jmk, i, j, l, nseg;
    int mj, im1, is1, is2, np1, lpp, twoK, threeK, nrwm1, wid;
    double vlen, wbrk[2], tnow;
    double x, y,w, rleny, temp, tc, dt, xm[2], ym[2], flm, deri2;
	double ts,te, machine_zero;

// ***** INITIAL SET. 
    // Parameter adjustments 
	twoK=k*2;
	threeK=k*3;

    // Function Body 
    klin = klini;
    if (klin != 1 && klin != 2)	klin=0;
	if(klin==0) wid=2; else wid=1;//id of weight coefficients.
    rleny = (float)npy/wind[3];
	rleny*=1.5;//increase resolutin by 50% since this is rational.
    np1 = n + 1;
    clip = wind[2] > 0.f;

    if(clip){
		if(nrw <= 1)
			return;
		x = wind[2] * .5f;
		xm[0] = wind[0] - x; xm[1] = wind[0] + x;
		y = wind[3] * .5f;
		ym[0] = wind[1] - y; ym[1] = wind[1] + y;
	}else{
		nrw = 2;
		rw[0] = t[k-1]; rw[1] = t[n];
	}

	machine_zero=bzmzro_();
// ******* CONVERT TO PP AND DRAW EACH LINE ***** 
    nrwm1=nrw-1;
    for (i = 1; i <= nrwm1; ++i){

	im1 = i-1;
	if(clip){
	    y=tc = (rw[im1] + rw[i]) * .5f;
//   *** TEST IF CURRENT LINE IN FRAME OR NOT 
	    w = bler_(k, n, t, rcoef[wid], tc, 0);
	    x = bler_(k, n, t, rcoef[0], tc, 0)/w;
	    if (klin == 0) y = bler_(k, n, t, rcoef[1], tc, 0)/w;
		else if (klin == 1) {y = x; x = tc;}
	    if(x<=xm[0] || x>=xm[1] || y<=ym[0] || y>=ym[1]) continue;
	}
// ***** DRAW LINE FROM RW(I) TO RW(I+1) 
	is1=bk2fli_(np1, t, rw[im1]);
	is2=bk2fli_(np1, t, rw[i]);	while(t[is2-1]==rw[i]) is2--;
	flm = rleny*(rw[i]-rw[im1])/(float)(is2-is1+1);
//   *** MOVE TO THE FIRST POSITION 
	y= tnow = rw[im1];
	w = bler_(k, n, t, rcoef[wid], tnow, 0);
	x = bler_(k, n, t, rcoef[0], tnow, 0)/w;
	if (klin == 0) y=bler_(k, n, t, rcoef[1], tnow, 0)/w;
	else if (klin == 1){y = x; x = tnow;}
	if(kfunc==1) (*movea2)((int)(x+.5),(int)(y+.5));
	else if(kfunc==2) (*movea2)((float)x,(float)y);
	else (*movea2)(x,y);

//   *** DRAW LINE ONE KNOT SPAN BY ONE 
	for (j = is1; j <= is2; ++j) {
		ts=t[j-1]; te=t[j];
		if(ts>=te) continue;
//     ...CONVERT THE ONE SPAN TO PP-REP 
	    jmk=j-k;
	    blcbp_(k,k, &t[jmk], &rcoef[0][jmk], &wk1[threeK], wbrk, wk1, &lpp);
	    blcbp_(k,k, &t[jmk], &rcoef[wid][jmk], &wk1[threeK], wbrk, wk1+twoK, &lpp);
	    w=wk1[twoK];
	    vlen=(wk1[1]-wk1[twoK+1]*wk1[0]/w)/w; if(vlen<0.) vlen=-vlen;
			//length of 1st derivative(about x coordinate).
	    if(klin == 0){
			blcbp_(k,k,&t[jmk],&rcoef[1][jmk],&wk1[threeK],wbrk,&wk1[k],&lpp);
			deri2=(wk1[k+1]-wk1[twoK+1]*wk1[k]/w)/w; if(deri2<0.) deri2=-deri2;
				//length of 1st derivative of y coordinate.
			if(vlen<deri2) vlen=deri2+vlen*.5;
			else           vlen=vlen+deri2*.5;
			//vlen is a good approximation of sqrt(vlen*vlen+deri2*deri2).
	    }//     ...NOW PP-REP IN WK1(.,1-2) 
// Computing MIN 
	    mj = (int)(flm*vlen) + 1;
		dt=(te-ts)/(float)mj;

		if(j==is1 || j==is2){
		    temp=rw[i]; if(temp>te) temp=te;
			if(ts>tnow) tnow=ts;
			if(dt<=machine_zero) nseg=1; else nseg=(int)((temp-tnow)/dt)+1;
			dt=(temp-tnow)/(float)nseg;

			if(j>is1){//Case that j is the last knot span and not the first one.
				//move to start point of the knot.
				y=tnow=ts; w=wk1[twoK]; x=wk1[0]/w;
				if(klin==0) y=wk1[k]/w;else if(klin==1){y=x; x=tnow;}
				if(kfunc==1) (*linea2)((int)(x+.5),(int)(y+.5));
				else if(kfunc==2) (*linea2)((float)x,(float)y);
				else (*linea2)(x,y);
			}
			//Draw middle points between the knots.
			for(l=1; l<nseg; l++){
				y=tnow+=dt;
				w=bpval2_(wbrk,wk1+twoK,k,tnow,1,0);
				x=bpval2_(wbrk,wk1,k,tnow,1,0)/w;
				if(klin==0) y=bpval2_(wbrk,wk1+k,k,tnow,1,0)/w;
				else if (klin == 1){y=x; x=tnow;}
				if(kfunc==1) (*linea2)((int)(x+.5),(int)(y+.5));
				else if(kfunc==2) (*linea2)((float)x,(float)y);
				else (*linea2)(x,y);
			}
		}else{
			//move to the start point of the knot.
			y=tnow=ts; w=wk1[twoK]; x=wk1[0]/w;
			if(klin==0) y=wk1[k]/w;else if(klin==1){y=x; x=tnow;}
			if(kfunc==1) (*linea2)((int)(x+.5),(int)(y+.5));
			else if(kfunc==2) (*linea2)((float)x,(float)y);
			else (*linea2)(x,y);

			//Draw middle points of between the knot BY INCREMENTING DT
			for(l=1; l<mj; l++){
				y=tnow+=dt;
				w=bpval2_(wbrk,wk1+twoK,k,tnow,1,0);
				x=bpval2_(wbrk,wk1,k,tnow,1,0)/w;
				if(klin==0) y=bpval2_(wbrk,&wk1[k],k,tnow,1,0)/w;
				else if (klin == 1){y=x; x=tnow;}
				if(kfunc==1) (*linea2)((int)(x+.5),(int)(y+.5));
				else if(kfunc==2) (*linea2)((float)x,(float)y);
				else (*linea2)(x,y);
			}
		}
	}

// DRAW THE LAST 1-SPAN. 
	y=tnow=rw[i];
	w=bpval2_(wbrk, wk1+twoK, k, tnow, 1, 0);
	x=bpval2_(wbrk, wk1, k, tnow, 1, 0)/w;
	if(klin==0) y=bpval2_(wbrk,&wk1[k],k,tnow,1,0)/w;
	else if(klin==1){ y=x; x=tnow;}
	if(kfunc==1) (*linea2)((int)(x+.5),(int)(y+.5));
	else if(kfunc==2) (*linea2)((float)x,(float)y);
	else (*linea2)(x,y);

    }
}
