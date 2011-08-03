/********************************************************************/
/* Copyright (c) 2011 DG Technologies Inc. and Yuzi Mizuno          */
/* All rights reserved.                                             */
/* Granted under the MIT license (see mg/MGCL.h for detail)         */
/********************************************************************/
//TO BE CALLED IN B L G S M T 
//FROM A PRACTICAL GUIDE TO SPLINE AND UPDATED BY Y.MIZUNO,9/7,'83 
//
//PUT DELX = X(.+1) - X(.)  INTO  V(.,4), 
//PUT THE THREE BANDS OF  Q-TRANSP*D  INTO  V(.,1-3), AND 
//PUT THE THREE BANDS OF (D*Q)-TRANSP*(D*Q)  AT AND ABOVE THE DIAGONAL 
//INTO  V(.,5-7) . 
//HERE,  Q IS  THE TRIDIAGONAL MATRIX OF ORDER (NPOINT-2,NPOINT) 
//WITH GENERAL ROW  1/DELX(I) , -1/DELX(I) - 1/DELX(I+1) , 1/DELX(I+1) 
//AND   D  IS THE DIAGONAL MATRIX  WITH GENERAL ROW  DY(I) . 
// *** INPUT *** 
// NCD..............SPACE DIMENSION.
// X[NPOINT]........ordinates of the data points.
// DY(I)............ERROR ESTIMATE, 1<=I<=NPOINT. 
// Y(NPOINT,J)......ORRESPONDING DATA ORDINATES at X[i] OF NCD SPACE DIMENSION
//                    (1<=J<=NCD). 
// NPOINT...........NUM OF POINTS. 
// *** OUTPUT *** 
// V(NPOINT,7) 
// QTY(NPOINT,NCD)..Q-TRANSP. * Y
void blgsm1_(int ncd, const double *x, const double *dy, 
	const double *y, int npoint, int iy, double *v,double *qty
);
