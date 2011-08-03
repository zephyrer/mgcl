/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
typedef int (*S_fp)();
#define min(a,b)            (((a) <= (b)) ? (a) : (b))
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
	const double *rcoefy, int nwk2, double *wk1, double *wk2, double *rw);
