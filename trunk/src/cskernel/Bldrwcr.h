/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
typedef int (*S_fp)();
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
			const double *t,const double **rcoef, double *wk1);
