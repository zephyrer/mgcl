/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
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
	double *w1, double *w2, int *nrw, double *rw);
