/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
// TO BE CALLED IN  B L G S M T 
// FROM A PRACTICAL GUIDE TO SPLINE BY C.DE-BOOR, AND 
// UPDATED BY Y. MIZUNO. 
//
// CONSTRUCTS THE UPPER THREE DIAGS. IN V(I,J), I=2,,NPOINT-1, J=1,3, OF 
//  THE MATRIX  6*(1-P)*Q-TRANSP.*(D**2)*Q + P*R, THEN COMPUTES ITS 
//  L*L-TRANSP. DECOMPOSITION AND STORES IT ALSO IN V, THEN APPLIES 
//  FORWARD AND BACK SUBSTITUTION TO THE RIGHT SIDE Q-TRANSP.*Y IN  QTY 
//  TO OBTAIN THE SOLUTION IN  U . 
// *** INPUT *** 
// P..................PARAMETER P. 
// V(NPOINT,7)........OUTPUT OF BLGSM1. 
// QTY(NPOINT,NCD)....Q-TRANSP. * Y (OUTPUT OF BLGSM1) 
// NPOINT.............NUM OF POINTS. 
// DY(I).........ERROR ESTIMATE, 1<=I<=NPOINT. 
// NCD................SPACE DIMENSION 
// *** OUTPUT *** 
// U(NPOINT,NCD) 
// QU(IQU,NCD)........Q*U 
// DQU................DY*QU. 
void blgsm2_(double p, double *v, const double *qty, 
	int npoint,const double *dy, int ncd, int iqu, 
	double *u, double *qu, double *dqu);
